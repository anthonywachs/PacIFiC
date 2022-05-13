#include "Grains_Exec.hh"
#include "MPIWrapperGrains.hh"
#include <sys/types.h>
#include <unistd.h>


MPIWrapperGrains* Grains_Exec::wrapper = NULL;
list<App*> Grains_Exec::m_allApp;
bool Grains_Exec::m_exception_Contact = false; 
bool Grains_Exec::m_exception_Deplacement = false;
bool Grains_Exec::m_exception_Simulation = false;
bool Grains_Exec::m_MPI = false;
string Grains_Exec::m_TIScheme = "SecondOrderLeapFrog";
bool Grains_Exec::m_SphereAsPolyParaview = false;
int Grains_Exec::m_MPI_verbose = 2;
string Grains_Exec::m_ReloadType = "new" ;
Vecteur Grains_Exec::m_vgravite = VecteurNul;
Vecteur* Grains_Exec::m_translationParaviewPostProcessing = NULL ;
bool Grains_Exec::m_withdemcfd = false;
// bool Grains_Exec::m_slipveloutput = false;
bool Grains_Exec::m_withlubrication = false;
bool Grains_Exec::m_ContactDissipation = false;
bool Grains_Exec::m_ContactforceOutput = false;
bool Grains_Exec::m_ContactforceOutput_instantaneous = false;
bool Grains_Exec::m_withLiftForce = false;
bool Grains_Exec::m_withHydroForce = false;
bool Grains_Exec::m_withFluidTemperature = false;
bool Grains_Exec::m_withSolidTemperature = false;
bool Grains_Exec::m_isGrainsCompFeatures = false;
bool Grains_Exec::m_isGrainsPorosity = false;
bool Grains_Exec::m_stressTensor = false;
bool Grains_Exec::m_particleStressTensor = false;
bool Grains_Exec::m_withPressureGradient = false;
bool Grains_Exec::m_withStochasticDrag = false;
bool Grains_Exec::m_withStochasticNusselt = false;
bool Grains_Exec::m_addedmass_demcfd = false;
double Grains_Exec::m_Prandtl = 0.;
string* Grains_Exec::debugString = NULL;
size_t Grains_Exec::m_nb_total_particules = 0;
list< pair<Point*,VertexBase *> > Grains_Exec::m_allPolytopeRefPointBase;
list<IndexArray*> Grains_Exec::m_allPolytopeNodeNeighbors;
list<IndexArray*> Grains_Exec::m_allPolytopeNodesIndex;
list<vector< vector<int> >*> Grains_Exec::m_allPolyhedronFacesConnectivity;

bool Grains_Exec::m_periodique = false;
bool Grains_Exec::m_MPIperiodique = false;
string Grains_Exec::m_ReloadDirectory = "";
string Grains_Exec::m_SaveDirectory = "";
set<string> Grains_Exec::m_additionalDataFiles;
bool Grains_Exec::m_writingModeHybrid = false;
string Grains_Exec::m_GRAINS_HOME = ".";
string Grains_Exec::m_reloadFile_suffix = "B";
bool Grains_Exec::m_withCohesion = false;
vector<Scalar> Grains_Exec::m_stressTensorDomain = vector<Scalar>( 6, 0. );


// ----------------------------------------------------------------------------
// Constructeur 
Grains_Exec::Grains_Exec()
{}




// ----------------------------------------------------------------------------
// Destructeur
Grains_Exec::~Grains_Exec()
{}




// ----------------------------------------------------------------------------
// Destruction de la memoire (garbage collector)
void Grains_Exec::GarbageCollector()
{
  if ( m_translationParaviewPostProcessing ) 
    delete m_translationParaviewPostProcessing;

  if ( !m_allPolytopeRefPointBase.empty() )
  {
    for (list< pair<Point*,VertexBase *> >::iterator 
    	ilPointBase=m_allPolytopeRefPointBase.begin(); 
  	ilPointBase!=m_allPolytopeRefPointBase.end(); ilPointBase++)
    {	  
      delete [] ilPointBase->first;
      delete ilPointBase->second;
    } 
    m_allPolytopeRefPointBase.clear();
  }
  
  if ( !m_allPolytopeNodeNeighbors.empty() )
  {
    for (list<IndexArray*>::iterator ilIA=m_allPolytopeNodeNeighbors.begin();
    	ilIA!=m_allPolytopeNodeNeighbors.end(); ilIA++)	  
      // les pointeurs pointent sur des tableaux d'IndexArray, d'ou
      // l'utilisation de "delete []" au lieu de "delete"
      delete [] *ilIA;
    m_allPolytopeNodeNeighbors.clear();
  }
  
  if ( !m_allPolytopeNodesIndex.empty() )
  {
    for (list<IndexArray*>::iterator ilNI=m_allPolytopeNodesIndex.begin();
    	ilNI!=m_allPolytopeNodesIndex.end(); ilNI++)	  
      delete *ilNI;
    m_allPolytopeNodeNeighbors.clear();
  } 
  
  if ( !m_allPolyhedronFacesConnectivity.empty() )
  {
    for (list<vector< vector<int> >*>::iterator 
    	ilFC=m_allPolyhedronFacesConnectivity.begin();
    	ilFC!=m_allPolyhedronFacesConnectivity.end(); ilFC++)	  
      delete *ilFC;
    m_allPolytopeNodeNeighbors.clear();
  }            
}




// ----------------------------------------------------------------------------
// Acces au wrapper MPI de Grains
MPIWrapperGrains* Grains_Exec::getComm()
{
  return wrapper;
}




// ----------------------------------------------------------------------------
// Wrapper MPI de Grains
void Grains_Exec::setComm( MPIWrapperGrains* wrapper_ )
{
  wrapper = wrapper_;
}




// ----------------------------------------------------------------------------
// Accessor to the list of applications
list<App*> Grains_Exec::get_listApp()
{
  return m_allApp;
}




// ----------------------------------------------------------------------------
// Accessor to the list of applications
void Grains_Exec::set_listApp( list<App*> allApp_ )
{
  m_allApp = allApp_;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture d'un double dans un format donné (nombre de caracteres) 
string Grains_Exec::doubleToString( const double &figure, const int &size )
{
  ostringstream oss;
  oss.width(size);
  oss << left << figure; 
  
  return ( oss.str() ); 
} 




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture d'un double dans un format donné (precision) 
string Grains_Exec::doubleToString( ios_base::fmtflags format, int digits,
      	double const& number )
{ 
  ostringstream oss;
  if ( number != 0. )
  {
    oss.setf( format, ios::floatfield );
    oss.precision( digits );
  }
  oss << number; 
   
  return ( oss.str() ); 

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture d'un integer dans un string
string Grains_Exec::intToString( const int &figure )
{
  ostringstream oss;
  oss << figure; 
  
  return ( oss.str() ); 
} 




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Memoire utilisee par la simulation sur ce processeur 
size_t Grains_Exec::used_memory( void )
{
  ostringstream os ;      
  string word ;
  size_t result = 0 ; 
   
  os << "/proc/" << getpid() << "/status" ;
      
  ifstream in( os.str().c_str() ) ;
  if( !in )
    cout << "Grains_Exec::used_memory : unable to open " << os.str() << endl ;
  else
  {
    while( !in.eof() )
    {
      in >> word ;
      if( word == "VmSize:" )
      {
        in >> result ;
        in >> word ;
        if( !( word == "kB" ) )
          cout << "Grains_Exec::used_memory : Unit is " << word << endl ;
        result *= 1000 ;
        break ;
      }
    }
    in.close() ;
  }

  return ( result ) ;
}



	
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la memoire utilisee par la simulation sur ce processeur
void Grains_Exec::display_memory( ostream& os, size_t memory )
{
  static size_t const mo = 1024*1024 ;
  static size_t const go = 1024*1024*1024 ;

  if( memory > go )
    os << ( (double) memory )/go << " Go" ;
  else if( memory > mo )
    os << ( (double) memory )/mo << " Mo" ;
  else os << memory << " octets" ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute une liste de points de reference utilises par les polytopes
void Grains_Exec::addOnePolytopeRefPointBase( Point* refPB, VertexBase *refVB )
{
  pair<Point*,VertexBase *> pp( refPB, refVB );
  m_allPolytopeRefPointBase.push_back( pp );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute une description des points sommets pour un type de polytope
void Grains_Exec::addOnePolytopeNodeNeighbors( IndexArray* idar )
{
  m_allPolytopeNodeNeighbors.push_back( idar );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute un tableau d'indice des sommets pour un type de polytope
void Grains_Exec::addOnePolytopeNodeIndex( IndexArray* idar )
{
  m_allPolytopeNodesIndex.push_back( idar );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute une connectivite des faces pour un type de polyedre
void Grains_Exec::addOnePolyhedronFaceConnectivity( 
	vector< vector<int> >* faceCon )
{
  m_allPolyhedronFacesConnectivity.push_back( faceCon );
}
 



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie le temps de la derniere ecriture dans un fichier
void Grains_Exec::checkTime_outputFile( const string& filename, 
    	double const& current_time )
{
  string tline, syscom;    
  istringstream iss;
  double last_output_time, tt;

  // Get last output time  
  ifstream fileIN( filename.c_str(), ios::in );
  if ( fileIN.is_open() )
  {
    getline( fileIN, tline );
    while ( !fileIN.eof() ) 
    { 
      iss.str( tline );
      iss >> last_output_time;
      iss.clear();
      getline( fileIN, tline );    
    }
    fileIN.close();
   
    // If last output time is greater than current time, delete all data
    // from last output time to current time 
    if ( last_output_time - current_time > 1.e-12 )
    {
      cout << "Time inconsistency in output file " 
    	<< filename << endl;
      ifstream fileIN_( filename.c_str(), ios::in ); 
      ofstream fileOUT( (filename+".tmp").c_str(), ios::out );
      getline( fileIN_, tline );
      while ( !fileIN_.eof() ) 
      { 
        iss.str( tline );
        iss >> tt;
        iss.clear();
        if ( current_time - tt > -1.e-10 )
          fileOUT << tline << endl;
        getline( fileIN_, tline );    
      }
      fileIN_.close();
      fileOUT.close(); 
    
      syscom = "mv " + filename + ".tmp" + " " + filename ;
      system(syscom.c_str());
    }      
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recupere le chemin a partir d'un nom complet de fichier
string Grains_Exec::extractRoot( string const& FileName )
{
  size_t pos = FileName.rfind( "/" );
  string root ;
  if ( pos != string::npos )
  { 
    root = FileName ;
    root.erase( root.begin() + pos, root.end() );
  }
  else root = "." ;
  
  return root ;
}  




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Extrait le nom du fichier uniquement a partir d'un nom complet de fichier
string Grains_Exec::extractFileName( string const& FileName )
{
  size_t pos = FileName.rfind( "/" );
  string result = FileName;
  if ( pos != string::npos )
    result.erase( result.begin(), result.begin() + pos + 1 );
  
  return result ;
}  




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie que tous les fichiers de reload sont bien dans le meme repertoire 
// de reload (fichiers polyedres et polygones essentiellement)
void Grains_Exec::checkAllFilesForReload()
{
  if ( !m_additionalDataFiles.empty() )
  {
    set<string>::iterator is;
    string fileName ;
    string cmd = "bash " + Grains_Exec::m_GRAINS_HOME 
    	+ "/ExecScripts/addFiles.exec";    
      
    for (is=m_additionalDataFiles.begin();is!=m_additionalDataFiles.end();is++)
    {
      fileName = Grains_Exec::extractFileName( *is );
      cmd += " " + *is + " " + fileName + " " + m_SaveDirectory; 
    } 

    system( cmd.c_str() );       
  }  
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recupere le nom du fichier de restart a partir du fichier RFTable.txt
string Grains_Exec::restartFileName_AorB( string const& rootName, 
  	string const& RFTable_ext )
{
  string buf, rfn;
  ifstream FILE_IN( ( rootName + RFTable_ext ).c_str(), ios::in );
  if ( FILE_IN.is_open() )
  {   
    while( !FILE_IN.eof() )
      FILE_IN >> buf >> rfn ;
  }
  else rfn = rootName;   
  FILE_IN.close() ; 
     
  return rfn;  
}
