#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Matlab_PostProcessingWriter.hh"
#include "Particule.H"
#include "Obstacle.H"
#include "Box.H"
#include "Cellule.H"
#include "Vecteur.H"
#include <zlib.h>
using namespace solid;

static int sizeof_Float32 = 4 ;
static int sizeof_Int32 = 4 ;


/* Constructeur par defaut 
--------------------------*/
Matlab_PostProcessingWriter::Matlab_PostProcessingWriter( DOMNode* dn,
	int const& rank_, int const& nbranks_ ):
  PostProcessingWriter( dn, rank_,nbranks_ ),
  m_MatlabReloadNumber( 0 ),
  m_binary( false ),
  BUFFER( NULL ),
  ALLOCATED( 0 ),
  OFFSET( 0 )
{
  m_MatlabFilename = ReaderXML::getNodeAttr_String( dn, "Name" );
  m_MatlabFilename_root = ReaderXML::getNodeAttr_String( dn, "Root" );
  string sm_binary = ReaderXML::getNodeAttr_String( dn, "Mode" );
  if ( sm_binary == "binary" ) m_binary = true;

}





/* Constructeur avec arguments 
------------------------------*/
Matlab_PostProcessingWriter::Matlab_PostProcessingWriter(
	int const& rank_,
  	int const& nbranks_,
	const string &name_,
	const string &root_,
  	const bool &isBinary):
  PostProcessingWriter( rank_, nbranks_ ),
  m_MatlabFilename_root( root_ ),
  m_MatlabFilename( name_ ),
  m_MatlabReloadNumber( 0 ),
  m_binary( isBinary ),
  BUFFER( NULL ),
  ALLOCATED( 0 ),
  OFFSET( 0 )
{}	




/* Destructeur 
--------------*/
Matlab_PostProcessingWriter::~Matlab_PostProcessingWriter()
{
}




/* Initialisation du post processeur 
------------------------------------*/
void Matlab_PostProcessingWriter::PostProcessing_start( Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC,
	vector<Fenetre> const& insert_windows )
{
  ostringstream ossRN;
  ossRN << m_MatlabReloadNumber;
  
  if ( Grains_Exec::m_ReloadType == "new" ) 
  {
    clearResultFiles();
    if ( Grains_Exec::m_MPI )
      Grains_Exec::getComm()->MPI_Barrier_ActivProc();
    
    if ( m_rank == 0 )
    { 
      string file;
      string partFilename = m_MatlabFilename + "_append_" 
	+ ossRN.str() + ".bin" ;
      file = m_MatlabFilename_root + "/" + partFilename ;
      m_Matlab_Append.open( file.c_str(), ios::out );
      
    } 
    
    one_output( temps, dt, particules );
  }
  else
  {
    // reload NOT IMPLEMENTED
    cout<< "Matlab_PPW grains reload NOT IMPLEMENTED" << endl;
  }
  
}




/* Relit un fichier binaire dans le cas d'un restart dans la continuite 
   et transfère son contenu dans le flux correspondant 
-------------------------------------------------------------------*/
void Matlab_PostProcessingWriter::readbinFile( string const& filename, 
	ostringstream& ossflux )
{
  // A Faire
}




/* Recupere le dernier numero de cycle dans le cas d'un restart 
   dans la continuite 
---------------------------------------------------------------*/
int Matlab_PostProcessingWriter::getPreviousCycleNumber() const
{
  int cyleNumber = 0;
  string tline, previous_tline, buffer, part;

  // Lecture dans le fichier pvd d'obstacles
  ifstream fileIN( (m_MatlabFilename_root + "/" + m_MatlabFilename
       	+ "_Obstacles.pvd" ).c_str(), ios::in );
  while ( tline != "</Collection>" )
  {
    previous_tline = tline;
    getline( fileIN, tline );
  }
  
  // Manipulation de la derniere ligne pour en extraire le numero de cycle
  istringstream iss( previous_tline ); 
  iss >> buffer >> buffer >> buffer >> buffer >> part;
  size_t pos = part.find( "_T" );
  string sub = part.substr( pos );
  sub.erase( sub.begin(), sub.begin() + 2 );
  sub.erase( sub.end() - 7, sub.end() ); 
  istringstream issNum( sub );
  issNum >> cyleNumber; 

  return cyleNumber;
}




/* Ecriture d'evolution 
-----------------------*/
void Matlab_PostProcessingWriter::PostProcessing( Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC )
{
  one_output( temps, dt, particules );  
}




/* Clot les ecritures 
---------------------*/
void Matlab_PostProcessingWriter::PostProcessing_end()
{
  if (m_rank == 0)
    m_Matlab_Append.close();
}




/* Ecriture
-----------*/
void Matlab_PostProcessingWriter::one_output( Scalar const& temps,
  	Scalar const& dt, list<Particule*> const* particules)
{
   
  list<Particule*>::const_iterator particule;
  Point gc; 
  Quaternion const* qrot; 
  //Vecteur const* vqrot;
  Vecteur const* vectrans; 
  Vecteur const* vecrot;
  Vecteur const* vecForce;
  Vecteur const* vecMoment;
  Vecteur const* PPTranslation = 
  	Grains_Exec::m_translationParaviewPostProcessing ;

  int nbpts=0;
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE && (*particule)->getTag() != 2 )
      nbpts++;
  
  // Position Xc - Yc - Zc - Ux - Vy - Wz - Ox - Oy - Oz
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 22*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE 
		&& (*particule)->getTag() != 2 )
      {
	//temps 
	write_double_binary( temps ) ;
	
	//Id particule
	write_double_binary(float((*particule)->getID()));
	
        //centre de gravite
	gc = *(*particule)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
	for (int comp=0;comp<3;++comp)
	  write_double_binary( gc[comp] ) ;

	//rayon particule
	write_double_binary((*particule)->getRayon());

	//rotation
	qrot  = (*particule)->getRotation();
	//vqrot = qrot->getVecteur();
	for (int comp=0;comp<4;++comp)
	  write_double_binary( (*qrot)[comp] ) ;	
 
	//vitesse de translation
	vectrans = (*particule)->getVitesseTranslation();
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vectrans)[comp] ) ;	
 
	//vitesse de rotation
	vecrot = (*particule)->getVitesseRotation();
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vecrot)[comp] ) ;

	// Force (Pas encore initialisee)
	vecForce = (*particule)->getForce();
	for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vecForce)[comp] ) ;
	
	// Moment (Pas encore initialise)
	vecMoment = (*particule)->getMoment();
	for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vecMoment)[comp] ) ;

      }
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE 
		&& (*particule)->getTag() != 2 )
      {
	//temps 
	m_Matlab_Append << temps << " " ;  
	
	//Id particule
	m_Matlab_Append << (*particule)->getID() << " " ; 

	//centre de gravite
        gc = *(*particule)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
	for (int comp=0;comp<3;++comp)
	  m_Matlab_Append << gc[comp] << " " ;  
	
	//rayon particule
	m_Matlab_Append << (*particule)->getRayon() << " " ; 

	//rotation
	qrot  = (*particule)->getRotation();
	//vqrot = qrot->getVecteur();
	for (int comp=0;comp<4;++comp)
	  m_Matlab_Append << (*qrot)[comp] << " " ;  
	
	//vitesse de translation
	vectrans = (*particule)->getVitesseTranslation();
        for (int comp=0;comp<3;++comp)
	  m_Matlab_Append << (*vectrans)[comp] << " " ;  
	
	//vitesse de rotation
	vecrot = (*particule)->getVitesseRotation();
        for (int comp=0;comp<3;++comp)
	  m_Matlab_Append << (*vecrot)[comp] << " " ;  
	
	// Force
	vecForce = (*particule)->getForce();
	for (int comp=0;comp<3;++comp)
	  m_Matlab_Append << (*vecForce)[comp] << " " ; 
	
	// Moment
	vecMoment = (*particule)->getMoment();
	for (int comp=0;comp<3;++comp)
	  m_Matlab_Append << (*vecMoment)[comp] << " " ;  
      }
    m_Matlab_Append << endl;
  }
    
  if ( m_binary )
  {
    m_Matlab_Append.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
  }  
  
   
}




/* Efface les fichiers resultats
--------------------------------*/
void Matlab_PostProcessingWriter::clearResultFiles() const
{
  if ( m_rank == 0 ) 
  {
    string cmd = "bash " + Grains_Exec::m_GRAINS_HOME 
     	+ "/ExecScripts/Matlab_clear.exec "; // + parameters
    cout << "Warning: Matlab_clear.exec not defined in $GRAINS_HOME/ExecScripts" 
    	<< endl;	
//    system( cmd.c_str() );
  }   
}

	

//-----------------------------------------------------------------------
// Methodes copiees de Pelicans-2.2.4: classe Matlab_PostProcessingWriter
//-----------------------------------------------------------------------
void Matlab_PostProcessingWriter:: start_output_binary( int size, int number )
{
  int current_output_size = size*number ;
  int ncomp = current_output_size + (current_output_size+999)/1000 
  	+ 12 + sizeof_Int32 ;
  check_allocated_binary( ncomp ) ;
  //CURRENT_LENGTH = store_int_binary(current_output_size) ;
}




void Matlab_PostProcessingWriter:: write_double_binary( double val )  
{
  *((float*)&(BUFFER[OFFSET])) = (float)val ;
  OFFSET += sizeof_Float32  ;
}



int Matlab_PostProcessingWriter:: store_int_binary( int val )  
{
  int result = OFFSET ;
  *((int*)&(BUFFER[OFFSET])) = val ;
  OFFSET += sizeof_Int32  ;
  return result ;
}




void Matlab_PostProcessingWriter:: check_allocated_binary( int size )  
{
  if(OFFSET+size>=ALLOCATED) 
  {
    int new_size = max( 2*ALLOCATED, (int)1024 ) ;
    new_size = max( new_size, 2*(OFFSET+size) ) ;
    new_size = 4 * ( new_size/4 +1 ) ; // allignement sur 4 bytes
      
    char * new_buffer = new char [ new_size ] ;
    for( int i=0 ;i<OFFSET ;i++ ) new_buffer[i] = BUFFER[i] ;
    if( BUFFER!=0 ) delete [] BUFFER ;
    BUFFER = new_buffer ;
    ALLOCATED = new_size ;      
  }
}


