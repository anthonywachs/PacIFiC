#include "GrainsMPI.H"
#include "GrainsMPITest.hh"
#include "Grains_BuilderFactory.H"
#include "GrainsParameters.H"
#include "GrainsPorosity.H"
#include "GrainsGeomTest.H"
#include "GrainsCoupledWithFluidMPI.hh"
#include "GrainsCompFeatures.H"
#include <string>
using namespace std;

EAPPLI Grains_BuilderFactory::m_context = UNDEFINED;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de l'application GRAINS avec l'option demandee
Grains* Grains_BuilderFactory::create(DOMElement* root)
{
  // Preconditions
  assert(root != NULL);

  Grains* grains = NULL;

  // Construction de l'application
  string type   = ReaderXML::getNodeName( root );
  string option = ReaderXML::getNodeAttr_String( root, "Type" );

  m_context = DIM_2;
  if ( type == "Grains3D" || type == "GrainsGeomTest" ) m_context = DIM_3;
  
  if ( option == "Standard" ) grains = new Grains();
  else if ( option == "MPI" ) grains = new GrainsMPI();
  else if ( option == "Parameters" ) grains = new GrainsParameters(); 
  else if ( option == "MPITest" ) grains = new GrainsMPITest(); 
  else if ( option == "VolumeInertie" ) grains = new GrainsCompFeatures();  
  else if ( option == "Porosity" ) grains = new GrainsPorosity();  
    
  // Postconditions
  assert( grains != NULL );

  return grains;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de l'application GrainsCoupledWithFluid avec l'option demandee
GrainsCoupledWithFluid* Grains_BuilderFactory::createCoupledWithFluid(
	DOMElement* root, double rhoFluide, double grid_size )
{
  // Preconditions
  assert( root != NULL );

  GrainsCoupledWithFluid* grains = NULL;

  // Construction de l'application
  string type   = ReaderXML::getNodeName( root );
  string option = ReaderXML::getNodeAttr_String( root, "Type" );

  m_context = DIM_2;
  if ( type == "Grains3D" ) m_context = DIM_3;

  if ( option == "CoupledFluid" ) 
    grains = new GrainsCoupledWithFluid( rhoFluide, grid_size );
  else if ( option == "CoupledFluidMPI" ) 
    grains = new GrainsCoupledWithFluidMPI( rhoFluide, grid_size );   

  // Postconditions
  assert( grains != NULL );

  return grains;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Contexte de la simulation
EAPPLI Grains_BuilderFactory::getContext()
{
  return m_context;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute le chemin vers les fichiers dtd en fonction de la variable
// shell GRAINS_HOME et du mode d'execution (dimension, avec ou sans 
// fluide, ... ), cree un fichier temporaire contenant ce chemin et renvoie le
// nom de ce fichier, c'est ce dernier qui sera lu par Grains
string Grains_BuilderFactory::init( string const& filename, 
	const int& rank, const int& nprocs )
{
  if ( rank == 0 )
  {        
    // Recupere la variable GRAINS_HOME et la stocke
    char* grainshome = getenv( "GRAINS_HOME" );
    string str_grainshome( grainshome );
    Grains_Exec::m_GRAINS_HOME = str_grainshome; 

    // Creation du fichier temporaire 
    string tline, buffer, 
  	header1 = "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>", 
	header2 = "<!DOCTYPE Grains3D SYSTEM \"",
	option, dtd_file ;
    list<string> inputFile_linelist;
    int dimension = 0;
  
    header2 += Grains_Exec::m_GRAINS_HOME + "/Main/dtd/" ;
   
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
  // Dans le cas d'un couplage avec le fluide, Grains3D peut tourner en
  // sequentiel alors que PeliGRIFF tourne en parallel  dans ce cas PeliGRIFF 
  // appelle cette methode avec nprocs = 1, ainsi MPI_Barrier n'est pas execute
  // (ce qui evite un blocage MPI où le master attendrait la reponse des autres
  // procs qui eux n'executent pas la commande)
  if ( nprocs > 1 ) MPI_Barrier( MPI_COMM_WORLD );  
  
  return( filename + ".tmp" );
} 
