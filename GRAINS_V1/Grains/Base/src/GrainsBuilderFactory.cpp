#include "GrainsBuilderFactory.hh"
#include "GrainsParameters.hh"
#include "GrainsCRBFeatures.hh"
#include "GrainsTestDev.hh"
#include <string>
using namespace std;

EAPPLI GrainsBuilderFactory::m_context = UNDEFINED;


// ----------------------------------------------------------------------------
// Creates and returns a standard Grains application
Grains* GrainsBuilderFactory::create( DOMElement* root )
{
  // Preconditions
  assert( root != NULL );

  Grains* grains = NULL;

  string type   = ReaderXML::getNodeName( root );
  string option = ReaderXML::getNodeAttr_String( root, "Type" );

  m_context = DIM_2;
  if ( type == "Grains3D" ) m_context = DIM_3;  
  
  if ( option == "Standard" ) grains = new Grains();
  else if ( option == "Parameters" ) grains = new GrainsParameters();
  else if ( option == "CompositeRigidBody" ) grains = new GrainsCRBFeatures();
  else if ( option == "TestDev" ) grains = new GrainsTestDev();     
    
  // Postconditions
  assert( grains != NULL );

  return ( grains );
}




// ----------------------------------------------------------------------------
// Returns the space dimension of the simulation
EAPPLI GrainsBuilderFactory::getContext()
{
  return ( m_context );
}




// ----------------------------------------------------------------------------
// Adds the path to the dtd files using the GRAINS_HOME variable to
// a copy of the input file. Returns the name of this copy
string GrainsBuilderFactory::init( string const& filename, int const& rank, 
  	int const& nprocs )
{
  // Store the input file name
  GrainsExec::m_inputFile = filename;
  
  // Get the GRAINS_HOME from the shell
  char* grainshome = getenv( "GRAINS_HOME" );
  string str_grainshome( grainshome );
  GrainsExec::m_GRAINS_HOME = str_grainshome;   
  
  if ( rank == 0 )
  {        
    // Creates the copy file
    string tline, buffer, 
  	header1 = "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>", 
	header2 = "<!DOCTYPE Grains3D SYSTEM \"",
	option, dtd_file ;
    list<string> inputFile_linelist;
    int dimension = 0;
  
    header2 += GrainsExec::m_GRAINS_HOME + "/Main/dtd/" ;
   
    ifstream fileIN( filename.c_str(), ios::in );
    ofstream fileOUT( ( filename + ".tmp" ).c_str(), ios::out );   

    while ( !fileIN.eof() ) 
    { 
      buffer.clear();
      getline( fileIN, tline, '\n' );
      istringstream iss(tline);
      iss >> buffer;
      if ( buffer != "<?xml" && buffer != "<!DOCTYPE" )  
      {
        inputFile_linelist.push_back( tline );
        if ( buffer == "<Grains3D" || buffer == "<GrainsGeomTest" ) 
        {
          dimension = 3;
	  buffer.clear();
          iss >> buffer;
	  size_t pos = buffer.find( "\"" );
          string sub = buffer.substr( pos + 1 );
	  pos = sub.rfind( "\"" );
	  option = sub.erase( pos );
        }	
        else if ( buffer == "<Grains2D" ) 
        {
          dimension = 2;
	  buffer.clear();
	  iss >> buffer;
	  size_t pos = buffer.find( "\"" );
          string sub = buffer.substr( pos + 1 );
	  pos = sub.rfind( "\"" );
	  option = sub.erase( pos );	
        }
        else if ( buffer == "<Grains3D>" || buffer == "<GrainsGeomTest>" )
          dimension = 3;
        else if ( buffer == "<Grains2D>" )
          dimension = 2;	      		
      }	
    }
  
    if ( dimension == 2 )
    {
      if ( option == "CoupledFluid" || option == "CoupledFluidMPI" )
        header2 += "Grains2D_InFluid.dtd\">";
      else
        header2 += "Grains2D.dtd\">";  
    }
    else
    {
      if ( option == "CoupledFluid" || option == "CoupledFluidMPI" )
        header2 += "Grains3D_InFluid.dtd\">";
      else
        header2 += "Grains3D.dtd\">";   
    }
 
    fileOUT << header1 << endl;
    fileOUT << header2 << endl;  
    for (list<string>::iterator il=inputFile_linelist.begin();
      il!=inputFile_linelist.end();il++)
      fileOUT << *il << endl;
  
    fileIN.close();
    fileOUT.close();
  }

  if ( nprocs > 1 ) MPI_Barrier( MPI_COMM_WORLD );  
  
  return ( filename + ".tmp" );
} 
