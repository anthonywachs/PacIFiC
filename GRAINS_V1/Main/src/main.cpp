#include <mpi.h>
#include "Grains.hh"
#include "GrainsBuilderFactory.hh"
#include "ReaderXML.hh"
#include <string>
using namespace std;


int main( int argc, char *argv[] )
{
  // MPI initialization
  MPI_Init(&argc,&argv);
  int rankproc = 0, nprocs = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );


  // Input file
  string filename = argv[1], filename_exe;
  size_t error = 0;
  size_t pos = filename.find(".xml");
  if ( pos == string::npos )
  {
    cout << "ERROR : input file need the .xml extension" << endl;
    error = 1;
  }


  // Execute the Grains application
  if ( !error )
  {
    // Create a temporary input file with the proper XML header
    filename_exe = GrainsBuilderFactory::init( filename, rankproc, nprocs );

    // Creates the Grains application
    ReaderXML::initialize();
    DOMElement* rootNode = ReaderXML::getRoot( filename_exe );
    string option = ReaderXML::getNodeAttr_String( rootNode, "Type" );

    Grains* grains = NULL;
    grains = GrainsBuilderFactory::create( rootNode );

    // Initial output message
    grains->initialOutputMessage();

    // Tasks to perform before time-stepping
    grains->do_before_time_stepping( rootNode );
    ReaderXML::terminate();

    // Delete the temporary input file
    if ( rankproc == 0 )
    {
      string cmd = "/bin/rm " + filename_exe;
      GrainsExec::m_return_syscmd = system( cmd.c_str() );
    }

    // Run the simulation
    grains->Simulation();

    // Tasks to perform after time-stepping
    grains->do_after_time_stepping();

    // Delete the Grains application
    delete grains;
  }

  // Close all MPI apps
  MPI_Finalize();


  return(0);
}
