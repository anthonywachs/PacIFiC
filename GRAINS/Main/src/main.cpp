#include <mpi.h>
#include "Grains.H"
#include "Grains_BuilderFactory.H"

#include "ReaderXML.hh"

#include <string>
using namespace std;

void help()
{
  cout << "Using ! \n grains inputFile \n"
       << endl;
}

// ============================================================================

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  int rankproc = 0, nprocs = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  // Fichier de mise en donnees
  string filename;
  if ( argc < 2 )
  {
    cout << "Filename : " << flush;
    cin >> filename;
  }
  else filename = argv[1];

  if (filename == "-h")
  {
    help();
    exit(1);
  }

  // Construction de l'application
  string filename_exe ;
  size_t pos = filename.find(".xml");
  if ( pos == string::npos )
    cout << "ERROR : input file need the .xml extension" << endl;
  else
  {
    filename_exe = Grains_BuilderFactory::init( filename, rankproc, nprocs );

    ReaderXML::initialize();
    DOMElement* rootNode = ReaderXML::getRoot( filename_exe );
    string option = ReaderXML::getNodeAttr_String( rootNode, "Type" );


    Grains* grains = NULL;
    if ( option == "CoupledFluid" || option == "CoupledFluidMPI" )
      grains = Grains_BuilderFactory::createCoupledWithFluid( rootNode, 1000.,
      	0. );
    else
      grains = Grains_BuilderFactory::create( rootNode );
    grains->Construction( rootNode );
    grains->Forces( rootNode );
    grains->Chargement( rootNode );

    ReaderXML::terminate();
    if ( rankproc == 0 )
    {
      string cmd = "/bin/rm " + filename_exe;
      system( cmd.c_str() );
    }

    if ( option == "CoupledFluid" || option == "CoupledFluidMPI" )
      grains->InitialPostProcessing();

    grains->Simulation();

    if ( option == "CoupledFluid" || option == "CoupledFluidMPI" )
      grains->doPostProcessing();

    delete grains;
  }

  MPI_Finalize();

  return(0);
}
