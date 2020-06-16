#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "CompFeatures_PostProcessingWriter.hh"
#include "Particule.H"
#include "Obstacle.H"
#include "Box.H"
#include "Cellule.H"
#include "Vecteur.H"
#include <zlib.h>
#include <string>
#include <iostream>
using namespace solid;



/* Constructeur par defaut 
--------------------------*/
CompFeatures_PostProcessingWriter::CompFeatures_PostProcessingWriter( DOMNode* dn,
	int const& rank_, int const& nbranks_ ):
  PostProcessingWriter( dn, rank_,nbranks_ ),
  m_ParaviewCycleNumber( 0 ),
  m_binary( false ),
  m_postProcessObstacle( true ),
  BUFFER( NULL ),
  ALLOCATED( 0 ),
  OFFSET( 0 )
{
  m_ParaviewFilename = ReaderXML::getNodeAttr_String( dn, "Name" );
  m_ParaviewFilename_root = ReaderXML::getNodeAttr_String( dn, "Root" );
  m_ParaviewCycleNumber = ReaderXML::getNodeAttr_Int( dn, "InitialCycleNumber" );
  string sm_binary = ReaderXML::getNodeAttr_String( dn, "Mode" );
  if ( sm_binary == "binary" ) m_binary = true;
  if ( ReaderXML::hasNodeAttr_String( dn, "Obstacle" ) )
  { 
    string sm_obstacle = ReaderXML::getNodeAttr_String( dn, "Obstacle" );
    if ( sm_obstacle == "no" ) m_postProcessObstacle = false; 
  } 
}





/* Constructeur avec arguments 
------------------------------*/
CompFeatures_PostProcessingWriter::CompFeatures_PostProcessingWriter(
	int const& rank_,
  	int const& nbranks_,
	const string &name_,
	const string &root_,
  	const bool &isBinary):
  PostProcessingWriter( rank_, nbranks_ ),
  m_ParaviewFilename_root( root_ ),
  m_ParaviewFilename( name_ ),
  m_ParaviewCycleNumber( 0 ),
  m_binary( isBinary ),
  m_postProcessObstacle( true ),  
  BUFFER( NULL ),
  ALLOCATED( 0 ),
  OFFSET( 0 )
{}	




/* Destructeur 
--------------*/
CompFeatures_PostProcessingWriter::~CompFeatures_PostProcessingWriter()
{
  vector<ostringstream*>::iterator iv;
  for (iv=m_Paraview_saveParticules_pvd.begin(); 
  	iv!=m_Paraview_saveParticules_pvd.end();iv++) delete *iv;
  m_Paraview_saveParticules_pvd.clear(); 
}




/* Initialisation du post processeur 
------------------------------------*/
void CompFeatures_PostProcessingWriter::PostProcessing_start( Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC,
	vector<Fenetre> const& insert_windows )
{
  size_t nbParticulesClasses = ParticuleClassesReference->size();

    clearResultFiles();
    if ( Grains_Exec::m_MPI )
      Grains_Exec::getComm()->MPI_Barrier_ActivProc();
    
    if ( m_rank == 0 )
    {
      // Particules
      ostringstream *ossNULL = NULL;
      m_Paraview_saveParticules_pvd.reserve( nbParticulesClasses );
      for (size_t i=0;i<nbParticulesClasses;++i)
        m_Paraview_saveParticules_pvd.push_back(ossNULL);
      for (size_t i=0;i<nbParticulesClasses;++i)
        m_Paraview_saveParticules_pvd[i] = new ostringstream;      
    
      for (size_t i=0;i<nbParticulesClasses;++i)
      {
        *m_Paraview_saveParticules_pvd[i] << "<?xml version=\"1.0\"?>" << endl;
        *m_Paraview_saveParticules_pvd[i] << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
        
	if ( m_binary ) 
          *m_Paraview_saveParticules_pvd[i] 
		<< " compressor=\"vtkZLibDataCompressor\"";
			
        *m_Paraview_saveParticules_pvd[i] << ">" << endl;
        *m_Paraview_saveParticules_pvd[i] << "<Collection>" << endl;
      }
    
    } 

    one_output_features( pwait );
}




/* Relit un fichier pvd dans le cas d'un restart dans la continuite 
   et transfère son contenu dans le flux correspondant 
-------------------------------------------------------------------*/
void CompFeatures_PostProcessingWriter::readPVDFile( string const& filename, 
	ostringstream& ossflux )
{
}




/* Recupere le dernier numero de cycle dans le cas d'un restart 
   dans la continuite 
---------------------------------------------------------------*/
int CompFeatures_PostProcessingWriter::getPreviousCycleNumber() const
{
  return(0);
}




/* Ecriture d'evolution 
-----------------------*/
void CompFeatures_PostProcessingWriter::PostProcessing( Scalar const& temps, 
	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC )
{
  one_output_features( pwait );  
}




/* Ecriture
-----------*/
void CompFeatures_PostProcessingWriter::one_output_features( 
	list<Particule*> const* pwait )
{
  list<Particule*>::const_iterator particule;
  list<string> Scalars;
  Scalars.push_back("NormU");
  Scalars.push_back("NormOm");
  Scalars.push_back("CoordNumb");   
  ostringstream ossCN, ossRK;
  ossCN << m_ParaviewCycleNumber; 
  ossRK << m_rank;   
  Scalar temps = 0.;

  string partFilename = m_ParaviewFilename + "_Particles_T" + ossCN.str();
  *m_Paraview_saveParticules_pvd[0] << "<DataSet timestep=\"" << temps 
  	<< "\" " << "group=\"\" part=\"0\" file=\"" << partFilename 
   << ".pvtu\"/>" << endl;       
  
  ofstream g( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
  	+ "_Particles.pvd" ).c_str(), ios::out );
  g << m_Paraview_saveParticules_pvd[0]->str();
  g << "</Collection>" << endl;
  g << "</VTKFile>" << endl;
  g.close();

  writePVTU_Paraview( partFilename, &empty_string_list,
        	&empty_string_list, &Scalars );     

  writeParticulesPostProcessing_Paraview( pwait,
   	partFilename + "_" + ossRK.str() + ".vtu" );

  m_ParaviewCycleNumber++; 
}



/* Clot les ecritures 
---------------------*/
void CompFeatures_PostProcessingWriter::PostProcessing_end()
{}




/* Ecriture
-----------*/
void CompFeatures_PostProcessingWriter::one_output( Scalar const& temps,
  	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
  	Obstacle *obstacle,
	LinkedCell const* LC )
{
}




/* Ecriture des obstacles
-------------------------*/
void CompFeatures_PostProcessingWriter::writeObstaclesPostProcessing_Paraview(
	list<MonObstacle*> const &allObstacles, string const &obsFilename )
{
} 




/* Mise a jour de l'indicateur des obstacles 
--------------------------------------------*/
void CompFeatures_PostProcessingWriter::updateObstaclesIndicator(
	Scalar const& temps,
  	Scalar const& dt,
	Obstacle *obstacle )
{
}




/* Ecriture des particules
--------------------------*/
void CompFeatures_PostProcessingWriter::writeParticulesPostProcessing_Paraview(
	list<Particule*> const* particules, const string &partFilename,
	bool const& forceForAllTag )
{
  list<Particule*>::const_iterator particule;
  Vecteur const* PPTranslation = 
  	Grains_Exec::m_translationParaviewPostProcessing ;

  ofstream f( ( m_ParaviewFilename_root + "/" + partFilename ).c_str(), 
  	ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0, i;
  for (particule=particules->begin();particule!=particules->end();particule++)
  {
      nbpts += (*particule)->numberOfPoints_PARAVIEW();
      nbcells += (*particule)->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;

  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (particule=particules->begin();particule!=particules->end();particule++)
      (*particule)->write_polygonsPts_PARAVIEW( f, PPTranslation );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity, offsets, cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber = 0, last_offset = 0;
  for (particule=particules->begin();particule!=particules->end();particule++)
      (*particule)->write_polygonsStr_PARAVIEW(connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (ii=connectivity.begin();ii!=connectivity.end();ii++)
    f << *ii << " ";	
  f << endl;      
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (ii=offsets.begin();ii!=offsets.end();ii++)
    f << *ii << " ";	
  f << endl; 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (ii=cellstype.begin();ii!=cellstype.end();ii++)
    f << *ii << " ";	
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;

  f << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">" << endl;
  f << "<DataArray type=\"Float32\" Name=\"NormU\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (particule=particules->begin();particule!=particules->end();particule++)
    {
      double normU = Norm( *(*particule)->getVitesseTranslation() );
      int nc = (*particule)->numberOfCells_PARAVIEW();
      for (i=0;i<nc;++i) f << normU << " ";
    }
  f << endl;
  f << "</DataArray>" << endl;      

  f << "<DataArray type=\"Float32\" Name=\"NormOm\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (particule=particules->begin();particule!=particules->end();particule++)
    {
      double normOm = Norm( *(*particule)->getVitesseRotation() );
      int nc = (*particule)->numberOfCells_PARAVIEW();
      for (i=0;i<nc;++i) f << normOm << " ";
    }
  f << endl;
  f << "</DataArray>" << endl; 

  f << "<DataArray type=\"Float32\" Name=\"CoordNumb\" ";
  f << "format=\"ascii\">";
  f << endl;
  for (particule=particules->begin();particule!=particules->end();particule++)
    {
      double coordNum = double((*particule)->getCoordinationNumber());
      int nc = (*particule)->numberOfCells_PARAVIEW();
      for (i=0;i<nc;++i) f << coordNum << " ";
    }
  f << endl;
  f << "</DataArray>" << endl; 
  f << "</CellData>" << endl;                 
  
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  f << "</VTKFile>" << endl;	
  f.close();	    
}	




/* Ecriture du maillage de cellules du LinkedCell
-------------------------------------------------*/
void CompFeatures_PostProcessingWriter::writeLinkedCellPostProcessing_Paraview(
	LinkedCell const* LC, const string &partFilename )
{
}




/* Ecriture des fenetres d'insertion
------------------------------------*/
void CompFeatures_PostProcessingWriter::writeInsertionPostProcessing_Paraview(
   	vector<Fenetre> const& insert_windows,
  	const string &partFilename )
{
}




/* Ecriture des vecteurs vitesse de translation & rotation des particules
-------------------------------------------------------------------------*/
void CompFeatures_PostProcessingWriter::
	writeVectorsMotionPostProcessing_Paraview(
	list<Particule*> const* particules, const string &partFilename )
{
}	




/* Ecriture des vecteurs force de contact entre composants
----------------------------------------------------------*/
void CompFeatures_PostProcessingWriter::writeVectorsForcePostProcessing_Paraview(
	list<Particule*> const* particules,
  	LinkedCell const* LC, const string &partFilename )
{
}	




/* Ecriture des particules de forme sphérique sous forme d'un vecteur
---------------------------------------------------------------------*/
void CompFeatures_PostProcessingWriter::
	writeSpheresPostProcessing_Paraview(
	list<Particule*> const* particules, const string &partFilename,
	bool const& forceForAllTag )
{
}
	



/* Ecriture du fichier pvtu 
---------------------------*/
void CompFeatures_PostProcessingWriter::writePVTU_Paraview( const string &filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar )
{
  list<string>::const_iterator il;

  ofstream f( ( m_ParaviewFilename_root + "/" + filename + ".pvtu" ).c_str(),
  	ios::out );
  f << "<?xml version=\"1.0\"?>" << endl; 
  f << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;  
  f << "<PPoints>" << endl;
  f << "<PDataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">"
  	<< endl;
  f << "</PDataArray>" << endl;
  f << "</PPoints>" << endl;  
  f << "<PCells>" << endl;
  f << "<PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">"
  	<< endl;
  f << "</PDataArray>" << endl;
  f << "<PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">" << endl;
  f << "</PDataArray>" << endl;
  f << "<PDataArray Name=\"types\" type=\"Int32\" format=\"ascii\">" << endl;
  f << "</PDataArray>" << endl;
  f << "</PCells>" << endl;    
  if ( pointVector->size() || pointScalar->size())
  {
    f << "<PPointData";
    if ( pointVector->size() )
    {
      il = pointVector->begin();
      f << " Vectors=\"" << *il;
      il++;
      for ( ;il!=pointVector->end();il++) f << "," << *il;
      f << "\"";      
    }
    if ( pointScalar->size() )
    {
      il = pointScalar->begin();
      f << " Scalars=\"" << *il;
      il++;
      for ( ;il!=pointScalar->end();il++) f << "," << *il;
      f << "\"";      
    }
    f << ">" << endl;
    for (il = pointVector->begin();il!=pointVector->end();il++)
    {    
      f << "<PDataArray Name=\"" << *il 
    	<< "\" NumberOfComponents=\"3\" type=\"Float32\""
    	<< " format=\"ascii\">" << endl;
      f << "</PDataArray>" << endl;
    }
    for (il = pointScalar->begin();il!=pointScalar->end();il++)
    {    
      f << "<PDataArray Name=\"" << *il 
    	<< "\" type=\"Float32\" format=\"ascii\">" << endl;
      f << "</PDataArray>" << endl;
    }    
    f << "</PPointData>" << endl;    
  }      
  if ( cellScalar->size() )
  {        
    f << "<PCellData Scalars=\"";
    il = cellScalar->begin();
    f << *il;
    il++;
    for ( ;il!=cellScalar->end();il++) f << "," << *il;
    f << "\">"; 
    for (il = cellScalar->begin();il!=cellScalar->end();il++)
    {     
      f << "<PDataArray Name=\"" << *il << "\" type=\"Float32\"" 
    	<< " format=\"ascii\">" << endl;
      f << "</PDataArray>" << endl;
    }
    f << "</PCellData>" << endl;  
  }  
  for (int i=0;i<m_nprocs;++i)
  {
    f << "<Piece Source=\"" << filename << "_" << i << ".vtu\">" << endl;
    f << "</Piece>" << endl; 
  }
  f << "</PUnstructuredGrid>" << endl;
  f << "</VTKFile>" << endl;  
  f.close();  
}




/* Ecriture des composants pour une erreur de contact
-----------------------------------------------------*/
void CompFeatures_PostProcessingWriter::writeErreurComposantsPostProcessing(
	string const& filename,
  	list<Composant*> const& errcomposants )
{
}




/* Efface les fichiers resultats
--------------------------------*/
void CompFeatures_PostProcessingWriter::clearResultFiles() const
{
  ofstream clear_file("Paraview_clear_exec",ios::out);
  clear_file << "#!/bin/sh" << endl;
  clear_file << "# Paraview files" << endl;
  clear_file << "for file in " << m_ParaviewFilename_root << "/" <<
  	m_ParaviewFilename << "*.pvd" << endl;
  clear_file << "  do" << endl;
  clear_file << "    if [ -f $file ]" << endl;
  clear_file << "      then" << endl;
  clear_file << "        rm -f $file" << endl;
  clear_file << "      fi" << endl;
  clear_file << "    done" << endl;
  clear_file << endl;

  clear_file << "for file in " << m_ParaviewFilename_root << "/" <<
  	m_ParaviewFilename << "*.vtu" << endl;
  clear_file << "  do" << endl;
  clear_file << "    if [ -f $file ]" << endl;
  clear_file << "      then" << endl;
  clear_file << "        rm -f $file" << endl;
  clear_file << "      fi" << endl;
  clear_file << "    done" << endl;
  clear_file << endl;
  
  clear_file << "for file in " << m_ParaviewFilename_root << "/" <<
  	m_ParaviewFilename << "*.pvtu" << endl;
  clear_file << "  do" << endl;
  clear_file << "    if [ -f $file ]" << endl;
  clear_file << "      then" << endl;
  clear_file << "        rm -f $file" << endl;
  clear_file << "      fi" << endl;
  clear_file << "    done" << endl;
  clear_file << endl;    
  
  clear_file << "for file in " << m_ParaviewFilename_root << "/" <<
  	"ErreurContact*" << endl;
  clear_file << "  do" << endl;
  clear_file << "    if [ -f $file ]" << endl;
  clear_file << "      then" << endl;
  clear_file << "        rm -f $file" << endl;
  clear_file << "      fi" << endl;
  clear_file << "    done" << endl;
  clear_file << endl;  
  
  clear_file << "for file in " << m_ParaviewFilename_root << "/" <<
  	"LinkedCell*" << endl;
  clear_file << "  do" << endl;
  clear_file << "    if [ -f $file ]" << endl;
  clear_file << "      then" << endl;
  clear_file << "        rm -f $file" << endl;
  clear_file << "      fi" << endl;
  clear_file << "    done" << endl;
  clear_file << endl; 
  clear_file.close();

  system("chmod 740 Paraview_clear_exec");
  system("./Paraview_clear_exec"); 
  system("/bin/rm -f Paraview_clear_exec"); 
}

	



//-----------------------------------------------------------------------
// Methodes copiees de Pelicans-2.2.4: classe CompFeatures_PostProcessingWriter
//-----------------------------------------------------------------------
void CompFeatures_PostProcessingWriter:: start_output_binary( int size, int number )
{
}




void CompFeatures_PostProcessingWriter:: write_double_binary( double val )  
{
}




void CompFeatures_PostProcessingWriter:: write_int_binary( int val )  
{
}




int CompFeatures_PostProcessingWriter:: store_int_binary( int val )  
{
  return(0);
}




void CompFeatures_PostProcessingWriter:: check_allocated_binary( int size )  
{
}




void CompFeatures_PostProcessingWriter:: flush_binary( std::ofstream& file )  
{
}




void CompFeatures_PostProcessingWriter:: compress_segment_binary( int seg )  
{
}
