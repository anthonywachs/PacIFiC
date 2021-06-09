#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Paraview_PostProcessingWriter.hh"
#include "Particule.H"
#include "Composant.H"
#include "Obstacle.H"
#include "Box.H"
#include "Cylinder.H"
#include "Segment.H"
#include "Cellule.H"
#include "Vecteur.H"
#include <zlib.h>
using namespace solid;

static int sizeof_Float32 = 4 ;
static int sizeof_Int32 = 4 ;


/* Constructeur par defaut
--------------------------*/
Paraview_PostProcessingWriter::Paraview_PostProcessingWriter( DOMNode* dn,
    int const& rank_, int const& nbranks_ ):
  PostProcessingWriter( dn, rank_,nbranks_ ),
  m_ParaviewCycleNumber( 0 ),
  m_binary( false ),
  m_postProcessObstacle( true ),
  m_initialCycleNumber_forced( false ),
  m_network( false ),
  BUFFER( NULL ),
  ALLOCATED( 0 ),
  OFFSET( 0 )
{
  m_ParaviewFilename = ReaderXML::getNodeAttr_String( dn, "Name" );
  m_ParaviewFilename_root = ReaderXML::getNodeAttr_String( dn, "Root" );
  m_ParaviewCycleNumber = ReaderXML::getNodeAttr_Int( dn, "InitialCycleNumber" );
  string sm_binary = ReaderXML::getNodeAttr_String( dn, "Mode" );
  if ( ReaderXML::hasNodeAttr_String( dn, "Network" ) ) m_network = true;
  if ( sm_binary == "binary" ) m_binary = true;
  if ( ReaderXML::hasNodeAttr_String( dn, "Obstacle" ) )
  {
    string sm_obstacle = ReaderXML::getNodeAttr_String( dn, "Obstacle" );
    if ( sm_obstacle == "no" ) m_postProcessObstacle = false;
  }
}





/* Constructeur avec arguments
------------------------------*/
Paraview_PostProcessingWriter::Paraview_PostProcessingWriter(
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
  m_initialCycleNumber_forced( false ),
  m_network( false ),
  BUFFER( NULL ),
  ALLOCATED( 0 ),
  OFFSET( 0 )
{}




/* Destructeur
--------------*/
Paraview_PostProcessingWriter::~Paraview_PostProcessingWriter()
{
  vector<ostringstream*>::iterator iv;
  for (iv=m_Paraview_saveParticules_pvd.begin();
  	iv!=m_Paraview_saveParticules_pvd.end();iv++) delete *iv;
  m_Paraview_saveParticules_pvd.clear();
}




/* Initialisation du post processeur
------------------------------------*/
void Paraview_PostProcessingWriter::PostProcessing_start(
    Scalar const& temps,
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

  if ( Grains_Exec::m_ReloadType == "new" )
  {
    clearResultFiles();
    if ( Grains_Exec::m_MPI )
      Grains_Exec::getComm()->MPI_Barrier_ActivProc();

    if ( m_rank == 0 )
    {
      // Obstacles
      if ( m_postProcessObstacle )
      {
        m_Paraview_saveObstacles_pvd << "<?xml version=\"1.0\"?>" << endl;
        m_Paraview_saveObstacles_pvd <<
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
        if ( m_binary )
          m_Paraview_saveObstacles_pvd <<
	  	" compressor=\"vtkZLibDataCompressor\"";
        m_Paraview_saveObstacles_pvd << ">" << endl;
        m_Paraview_saveObstacles_pvd << "<Collection>" << endl;

        if ( Grains_Exec::m_periodique == true )
        {
          m_Paraview_saveObstaclesPeriodiques_pvd << "<?xml version=\"1.0\"?>"
      		<< endl;
          m_Paraview_saveObstaclesPeriodiques_pvd <<
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";

	  if ( m_binary )
            m_Paraview_saveObstaclesPeriodiques_pvd
      		<< " compressor=\"vtkZLibDataCompressor\"";

          m_Paraview_saveObstaclesPeriodiques_pvd << ">" << endl;
          m_Paraview_saveObstaclesPeriodiques_pvd << "<Collection>" << endl;
        }
      }

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

      if ( Grains_Exec::m_periodique || Grains_Exec::m_MPIperiodique )
      {
        m_Paraview_saveClonesPeriodiques_pvd << "<?xml version=\"1.0\"?>" << endl;
        m_Paraview_saveClonesPeriodiques_pvd <<
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";

	if ( m_binary )
          m_Paraview_saveClonesPeriodiques_pvd
		<< " compressor=\"vtkZLibDataCompressor\"";

        m_Paraview_saveClonesPeriodiques_pvd << ">" << endl;
        m_Paraview_saveClonesPeriodiques_pvd << "<Collection>" << endl;
      }

      if ( LC )
      {
        // Vitesse de translation
        m_Paraview_saveVectors_pvd << "<?xml version=\"1.0\"?>" << endl;
        m_Paraview_saveVectors_pvd <<
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";

        if ( m_binary )
          m_Paraview_saveVectors_pvd
		<< " compressor=\"vtkZLibDataCompressor\"";

        m_Paraview_saveVectors_pvd << ">" << endl;
        m_Paraview_saveVectors_pvd << "<Collection>" << endl;

        // Forces de contact
        m_Paraview_saveContactForces_pvd << "<?xml version=\"1.0\"?>" << endl;
        m_Paraview_saveContactForces_pvd <<
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";

        if ( m_binary )
          m_Paraview_saveContactForces_pvd
		<< " compressor=\"vtkZLibDataCompressor\"";

        m_Paraview_saveContactForces_pvd << ">" << endl;
        m_Paraview_saveContactForces_pvd << "<Collection>" << endl;

	if ( m_network )
	{
	  // Reseau Forces de contact
	  m_Paraview_saveForceChain_pvd << "<?xml version=\"1.0\"?>" << endl;
	  m_Paraview_saveForceChain_pvd <<
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";

	  if ( m_binary )
	    m_Paraview_saveForceChain_pvd
		<< " compressor=\"vtkZLibDataCompressor\"";

	  m_Paraview_saveForceChain_pvd<< ">" << endl;
	  m_Paraview_saveForceChain_pvd << "<Collection>" << endl;
	}
        // Maillage LinkedCell
        writePVTU_Paraview( "LinkedCell" , &empty_string_list,
		&empty_string_list, &empty_string_list );
      }
    }

     // Maillage LinkedCell
//      if ( LC )
//      {
//        ostringstream ossRK;
//        ossRK << m_rank;
//        writeLinkedCellPostProcessing_Paraview( LC,
//       	"LinkedCell_" + ossRK.str() + ".vtu" );
//      }

    one_output( temps, dt, particules, pperiodiques, ParticuleClassesReference,
  	obstacle, LC );
  }
  else
  {
    // En reload "same", l'existence du fichier de post-processing des obstacles
    // est obligatoire, donc on force m_postProcessObstacle = true
    m_postProcessObstacle = true ;
    if ( m_rank == 0 )
    {
      // Obstacles
      readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_Obstacles.pvd", m_Paraview_saveObstacles_pvd );

      if ( Grains_Exec::m_periodique == true )
        readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_ObstaclesPeriodiques.pvd",
	m_Paraview_saveObstaclesPeriodiques_pvd );

      // Particules
      ostringstream *ossNULL = NULL;
      m_Paraview_saveParticules_pvd.reserve( nbParticulesClasses );
      for (size_t i=0;i<nbParticulesClasses;++i)
        m_Paraview_saveParticules_pvd.push_back(ossNULL);
      for (size_t i=0;i<nbParticulesClasses;++i)
        m_Paraview_saveParticules_pvd[i] = new ostringstream;

      if ( nbParticulesClasses == 1 )
        readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_Particles.pvd", *m_Paraview_saveParticules_pvd[0] );
      else
        for (size_t i=0;i<nbParticulesClasses;++i)
	{
          ostringstream* ossPC = new ostringstream;
          *ossPC << i;
          readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
		+ "_Particles_Class"+ossPC->str() + ".pvd",
		*m_Paraview_saveParticules_pvd[i] );
	}

      if ( Grains_Exec::m_periodique || Grains_Exec::m_MPIperiodique )
        readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_ClonesPeriodiques.pvd", m_Paraview_saveClonesPeriodiques_pvd );

      if ( LC )
      {
        // Vitesse de translation
	readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_VecMotion.pvd", m_Paraview_saveVectors_pvd );

        // Forces de contact
	readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_VecContactForce.pvd", m_Paraview_saveContactForces_pvd );
      }

      if ( m_network )
      {
        // Reseau Forces de contact
	readPVDFile( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_ForceChain.pvd", m_Paraview_saveForceChain_pvd );
      }
    }

    // Cycle number
    if ( !m_initialCycleNumber_forced )
      m_ParaviewCycleNumber = getPreviousCycleNumber();
    ++m_ParaviewCycleNumber;
  }

  // Fenetres d'insertion
  if ( m_rank == 0 )
    writeInsertionPostProcessing_Paraview( insert_windows,
    	m_ParaviewFilename + "_InsertionWindows" );
}




/* Relit un fichier pvd dans le cas d'un restart dans la continuite
   et transf�re son contenu dans le flux correspondant
-------------------------------------------------------------------*/
void Paraview_PostProcessingWriter::readPVDFile( string const& filename,
	ostringstream& ossflux )
{
  string tline, buffer, part, keyword;
  bool keep = true ;
  int cyleNumber = 0;

  ifstream fileIN( filename.c_str(), ios::in );
  getline( fileIN,tline );
  while ( !fileIN.eof() )
  {
    if ( tline != "</Collection>" && tline != "</VTKFile>" )
    {
      keep = true ;
      if ( m_initialCycleNumber_forced )
      {
        istringstream iss( tline );
        iss >> keyword;
        if ( keyword == "<DataSet" )
        {
          iss >> buffer >> buffer >> buffer >> part;
          size_t pos = part.find( "_T" );
          string sub = part.substr( pos );
          sub.erase( sub.begin(), sub.begin() + 2 );
	  pos = sub.rfind( "." );
          sub.erase( sub.begin() + pos, sub.end() );
          istringstream issNum( sub );
          issNum >> cyleNumber;
          if ( cyleNumber > m_ParaviewCycleNumber ) keep = false ;
        }
      }
    }
    else keep = false ;
    if ( keep ) ossflux << tline << endl;
    getline( fileIN, tline );
  }
  fileIN.close();
}




/* Recupere le dernier numero de cycle dans le cas d'un restart
   dans la continuite
---------------------------------------------------------------*/
int Paraview_PostProcessingWriter::getPreviousCycleNumber() const
{
  int cyleNumber = 0;
  string tline, previous_tline, buffer, part;

  // Lecture dans le fichier pvd d'obstacles
  ifstream fileIN( (m_ParaviewFilename_root + "/" + m_ParaviewFilename
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
void Paraview_PostProcessingWriter::PostProcessing(
    Scalar const& temps,
    Scalar const& dt,
    list<Particule*> const* particules,
    list<Particule*> const* pwait,
    list<Particule*> const* pperiodiques,
    vector<Particule*> const* ParticuleClassesReference,
    Obstacle *obstacle,
    LinkedCell const* LC )
{
  one_output( temps, dt, particules, pperiodiques, ParticuleClassesReference,
  	obstacle, LC );
}




/* Clot les ecritures
---------------------*/
void Paraview_PostProcessingWriter::PostProcessing_end()
{}




/* Ecriture
-----------*/
void Paraview_PostProcessingWriter::one_output(
    Scalar const& temps,
    Scalar const& dt,
    list<Particule*> const* particules,
    list<Particule*> const* pperiodiques,
    vector<Particule*> const* ParticuleClassesReference,
    Obstacle *obstacle,
    LinkedCell const* LC )
{
  size_t nbParticulesClasses = ParticuleClassesReference->size();
  list<string> Scalars;
  Scalars.push_back("NormU");
  Scalars.push_back("NormOm");
  Scalars.push_back("CoordNumb");
  if( Grains_Exec::m_withFluidTemperature ||
      Grains_Exec::m_withSolidTemperature )
    Scalars.push_back("Temperature");
  ostringstream ossCN, ossRK;
  ossCN << m_ParaviewCycleNumber;
  ossRK << m_rank;
  // Obstacles
  if ( m_postProcessObstacle )
  {
    updateObstaclesIndicator( temps, dt, obstacle );
    if ( m_rank == 0 )
    {
      list<MonObstacle*> allObstacles = obstacle->getObstacles();
      list<MonObstacle*> allObstaclesPeriodiques;
      if ( Grains_Exec::m_periodique == true )
      {
        list<MonObstacle*>::iterator il;
        for (il=allObstacles.begin();il!=allObstacles.end();)
          if ( (*il)->materiau() == "periode" )
          {
            allObstaclesPeriodiques.push_back(*il);
            il = allObstacles.erase( il );
          }
          else il++;
      }

      string obsFilename = m_ParaviewFilename + "_Obstacles_T" + ossCN.str() +
   	".vtu";
      m_Paraview_saveObstacles_pvd << "<DataSet timestep=\"" << temps
      	<< "\" " << "group=\"\" part=\"0\" file=\"" << obsFilename << "\"/>\n";

      ofstream f( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_Obstacles.pvd" ).c_str(), ios::out );
//       f << m_Paraview_saveObstacles_pvd.str();
      writeBigOSS( f, m_Paraview_saveObstacles_pvd );
      f << "</Collection>" << endl;
      f << "</VTKFile>" << endl;
      f.close();

      if ( obstacle ) writeObstaclesPostProcessing_Paraview( allObstacles,
     		obsFilename );

      if ( Grains_Exec::m_periodique == true )
      {
        string obsPerFilename = m_ParaviewFilename + "_ObstaclesPeriodiques_T"
       		+ ossCN.str() + ".vtu";
        m_Paraview_saveObstaclesPeriodiques_pvd << "<DataSet timestep=\"" <<
       		temps << "\" " << "group=\"\" part=\"0\" file=\"" <<
		obsPerFilename << "\"/>" << endl;

        ofstream h( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       		+ "_ObstaclesPeriodiques.pvd" ).c_str(), ios::out );
//         h << m_Paraview_saveObstaclesPeriodiques_pvd.str();
        writeBigOSS( h, m_Paraview_saveObstaclesPeriodiques_pvd );
        h << "</Collection>" << endl;
        h << "</VTKFile>" << endl;
        h.close();

        if ( obstacle ) writeObstaclesPostProcessing_Paraview(
       		allObstaclesPeriodiques, obsPerFilename );
      }
    }
  }

  // Particules
  if ( nbParticulesClasses == 1 )
  {
    string partFilename = m_ParaviewFilename + "_Particles_T" + ossCN.str();
    if ( m_rank == 0 )
    {
      *m_Paraview_saveParticules_pvd[0] << "<DataSet timestep=\"" << temps
      	<< "\" " << "group=\"\" part=\"0\" file=\"" << partFilename
	<< ".pvtu\"/>" << endl;

      ofstream g( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_Particles.pvd" ).c_str(), ios::out );
//       g << m_Paraview_saveParticules_pvd[0]->str();
      writeBigOSS( g, *m_Paraview_saveParticules_pvd[0] );
      g << "</Collection>" << endl;
      g << "</VTKFile>" << endl;
      g.close();

      if ( (*ParticuleClassesReference)[0]->getForme()->getConvex()
       	->isSphere() && !Grains_Exec::m_SphereAsPolyParaview )
      {
        list<string> ptVec;
        ptVec.push_back("Orientation");
        writePVTU_Paraview( partFilename, &ptVec, &Scalars,
	 	&empty_string_list );
      }
      else writePVTU_Paraview( partFilename, &empty_string_list,
	 	&empty_string_list, &Scalars );
    }

//    if ( (*ParticuleClassesReference)[0]->getForme()->getConvex()
//	 == NULL )
//    {
//      writeParticulesPostProcessing_Paraview( particules,
//    	partFilename + "_" + ossRK.str() + ".vtu" );
//    }

    // Are the particles in the PostProcessingWindow ?
    if( PostProcessingWriter::m_bPPWindow[m_rank] )
    {
      if( (*ParticuleClassesReference)[0]->getForme()->getConvex()
          ->isSphere() && !Grains_Exec::m_SphereAsPolyParaview )
	writeSpheresPostProcessing_Paraview( particules,
     	  partFilename + "_" + ossRK.str() + ".vtu" );
      else
	writeParticulesPostProcessing_Paraview( particules,
    	  partFilename + "_" + ossRK.str() + ".vtu" );
    }

  }
  else
  {
    list<Particule*> empty_list;
    vector< list<Particule*> > partPerClasse( nbParticulesClasses,
     	empty_list );
    list<Particule*>::const_iterator particule;

    for (particule=particules->begin();particule!=particules->end();
     	particule++)
      partPerClasse[(*particule)->getParticuleClasse()].push_back(*particule);

    for (size_t i=0;i<nbParticulesClasses;++i)
    {
      ostringstream* ossPC = new ostringstream;
      *ossPC << i;
      string partFilename = m_ParaviewFilename + "_Particles_Class" +
      	ossPC->str() + "_T" + ossCN.str();
      if ( m_rank == 0 )
      {
        *m_Paraview_saveParticules_pvd[i] << "<DataSet timestep=\"" << temps
       		<< "\" " << "group=\"\" part=\"0\" file=\"" << partFilename
		<< ".pvtu\"/>" << endl;

        ofstream g( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       		+ "_Particles_Class" + ossPC->str() + ".pvd" ).c_str(),
		ios::out );
//         g << m_Paraview_saveParticules_pvd[i]->str();
	writeBigOSS( g, *m_Paraview_saveParticules_pvd[i] );
        g << "</Collection>" << endl;
        g << "</VTKFile>" << endl;
        g.close();

        if ( (*ParticuleClassesReference)[i]->getForme()->getConvex()
	 	->isSphere() && !Grains_Exec::m_SphereAsPolyParaview )
        {
          list<string> ptVec;
          ptVec.push_back("Orientation");
          writePVTU_Paraview( partFilename, &ptVec, &Scalars,
	   	&empty_string_list );
        }
        else writePVTU_Paraview( partFilename, &empty_string_list,
	 	&empty_string_list, &Scalars );
      }

      // Are the particles in the PostProcessingWindow ?
      if( PostProcessingWriter::m_bPPWindow[m_rank] )
      {
	if ( (*ParticuleClassesReference)[i]->getForme()->getConvex()
	 	  ->isSphere() && !Grains_Exec::m_SphereAsPolyParaview )
          writeSpheresPostProcessing_Paraview( &partPerClasse[i],
       		  partFilename + "_" + ossRK.str() + ".vtu" );
	else
          writeParticulesPostProcessing_Paraview( &partPerClasse[i],
       		  partFilename + "_" + ossRK.str() + ".vtu" );
	delete ossPC;
      }

    }
  }

  if ( Grains_Exec::m_periodique || Grains_Exec::m_MPIperiodique )
  {
    string partFilename = m_ParaviewFilename + "_ClonesPeriodiques_T" +
     	ossCN.str();
    if ( m_rank == 0 )
    {
      m_Paraview_saveClonesPeriodiques_pvd << "<DataSet timestep=\"" << temps
      	<< "\" " << "group=\"\" part=\"0\" file=\"" << partFilename
	<< ".pvtu\"/>" << endl;

      ofstream g( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_ClonesPeriodiques.pvd" ).c_str(), ios::out );
//       g << m_Paraview_saveClonesPeriodiques_pvd.str();
      writeBigOSS( g, m_Paraview_saveClonesPeriodiques_pvd );
      g << "</Collection>" << endl;
      g << "</VTKFile>" << endl;
      g.close();

      if ( nbParticulesClasses == 1
       	&& ( (*ParticuleClassesReference)[0]->getForme()->getConvex()
	->isSphere() && !Grains_Exec::m_SphereAsPolyParaview ) )
      {
        list<string> ptVec;
        ptVec.push_back("Orientation");
        writePVTU_Paraview( partFilename, &ptVec, &Scalars,
	 	&empty_string_list );
      }
      writePVTU_Paraview( partFilename, &empty_string_list, &empty_string_list,
       	&Scalars );
    }

    // Are the particles in the PostProcessingWindow ?
    if( PostProcessingWriter::m_bPPWindow[m_rank] )
    {
      if ( nbParticulesClasses == 1
      	  && ( (*ParticuleClassesReference)[0]->getForme()->getConvex()
	  ->isSphere() && !Grains_Exec::m_SphereAsPolyParaview ) )
	writeSpheresPostProcessing_Paraview( pperiodiques,
       		  partFilename + "_" + ossRK.str() + ".vtu",
		  Grains_Exec::m_MPIperiodique );
      else
	writeParticulesPostProcessing_Paraview( pperiodiques,
     		  partFilename + "_" + ossRK.str() + ".vtu",
		  Grains_Exec::m_MPIperiodique );
    }
  }

  if ( LC )
  {
    // Vitesse de translation & rotation
    string vectFilename = m_ParaviewFilename + "_VecMotion_T" + ossCN.str();
    if ( m_rank == 0 )
    {
      m_Paraview_saveVectors_pvd << "<DataSet timestep=\"" << temps
       	<< "\" " << "group=\"\" part=\"0\" file=\"" << vectFilename
	<< ".pvtu\"/>" << endl;

      ofstream h( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
       	+ "_VecMotion.pvd" ).c_str(), ios::out );
//       h << m_Paraview_saveVectors_pvd.str();
      writeBigOSS( h, m_Paraview_saveVectors_pvd );
      h << "</Collection>" << endl;
      h << "</VTKFile>" << endl;
      h.close();

      list<string> vecMotion;
      vecMotion.push_back("U");
      vecMotion.push_back("Omega");
      writePVTU_Paraview( vectFilename, &vecMotion, &empty_string_list,
       	&empty_string_list );
    }

    // Does this processor have to write down outputs ?
    if( PostProcessingWriter::m_bPPWindow[m_rank] )
      writeVectorsMotionPostProcessing_Paraview( particules,
   	vectFilename + "_" + ossRK.str() + ".vtu" );

    // Forces de contact
    string forceFilename = m_ParaviewFilename + "_VecContactForce_T" +
      	ossCN.str();
    if ( m_rank == 0 )
    {
      m_Paraview_saveContactForces_pvd << "<DataSet timestep=\"" << temps
        	<< "\" " << "group=\"\" part=\"0\" file=\"" << forceFilename
 	<< ".pvtu\"/>" << endl;

      ofstream h( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
        	+ "_VecContactForce.pvd" ).c_str(), ios::out );
//        h << m_Paraview_saveContactForces_pvd.str();
      writeBigOSS( h, m_Paraview_saveContactForces_pvd );
      h << "</Collection>" << endl;
      h << "</VTKFile>" << endl;
      h.close();

      list<string> vecForce;
      vecForce.push_back("Force");
      writePVTU_Paraview( forceFilename, &vecForce, &empty_string_list,
        	&empty_string_list );
    }

    // Does this processor have to write down outputs ?
    if( PostProcessingWriter::m_bPPWindow[m_rank] )
      writeVectorsForcePostProcessing_Paraview( particules, LC,
   	forceFilename + "_" + ossRK.str() + ".vtu", temps, dt );


    if ( m_network )
    {
      // Force chain
      string forceChainFilename = m_ParaviewFilename + "_ForceChain_T" +
        	ossCN.str();
      if ( m_rank == 0 )
      {
        m_Paraview_saveForceChain_pvd << "<DataSet timestep=\"" << temps
          	<< "\" " << "group=\"\" part=\"0\" file=\"" << forceChainFilename
          << ".pvtp\"/>" << endl;

        ofstream h( ( m_ParaviewFilename_root + "/" + m_ParaviewFilename
          	+ "_ForceChain.pvd" ).c_str(), ios::out );
        writeBigOSS( h, m_Paraview_saveForceChain_pvd );
        h << "</Collection>" << endl;
        h << "</VTKFile>" << endl;
        h.close();

        list<string> vecForce;
        vecForce.push_back("ForceChain");
        writePVTP_Paraview( forceChainFilename, &vecForce, &empty_string_list,
          	&empty_string_list );
      }

      // Does this processor have to write down outputs ?
      if( PostProcessingWriter::m_bPPWindow[m_rank] )
        writeForceChain_Paraview( particules, LC, forceChainFilename + "_" +
            ossRK.str() + ".vtp", temps, dt );
      }
  }

   m_ParaviewCycleNumber++;
}




/* Ecriture des obstacles
-------------------------*/
void Paraview_PostProcessingWriter::writeObstaclesPostProcessing_Paraview(
	list<MonObstacle*> const &allObstacles, string const &obsFilename )
{
  ofstream f( ( m_ParaviewFilename_root + "/" + obsFilename ).c_str(),
  	ios::out );
  list<MonObstacle*>::const_iterator il;

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0,i;
  for (il=allObstacles.begin();il!=allObstacles.end();il++)
  {
    //nbpts += (*il)->getForme()->getConvex()->numberOfPoints_PARAVIEW();
    //nbcells += (*il)->getForme()->getConvex()->numberOfCells_PARAVIEW();
    nbpts += (*il)->numberOfPoints_PARAVIEW();
    nbcells += (*il)->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point> ppp;
    list<Point>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (il=allObstacles.begin();il!=allObstacles.end();il++)
    {
      ppp = (*il)->getForme()->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/Points" );
  }
  else
    for (il=allObstacles.begin();il!=allObstacles.end();il++)
    {
      //(*il)->getForme()->write_polygonsPts_PARAVIEW( f );
      (*il)->write_polygonsPts_PARAVIEW( f );
    }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (il=allObstacles.begin();il!=allObstacles.end();il++)
    (*il)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );
  //  (*il)->getForme()->getConvex()->write_polygonsStr_PARAVIEW( connectivity,
  //  	offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/connectivity" );
  }
  else
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/offsets" );
  }
  else
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;

  f << "<CellData Scalars=\"Indicator\">" << endl;
  f << "<DataArray type=\"Float32\" Name=\"Indicator\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (il=allObstacles.begin();il!=allObstacles.end();il++)
  {
    double indic = (*il)->getIndicator();
    //int nc =(*il)->getForme()->getConvex()->numberOfCells_PARAVIEW();
    int nc =(*il)->numberOfCells_PARAVIEW();
    if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( indic );
    else for (i=0;i<nc;++i) f << indic << " ";
  }
  if ( m_binary ) flush_binary( f,
  	"writeObstaclesPostProcessing_Paraview/Indicator" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</CellData>" << endl;

  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();
}




/* Mise a jour de l'indicateur des obstacles
--------------------------------------------*/
void Paraview_PostProcessingWriter::updateObstaclesIndicator(
	Scalar const& temps,
  	Scalar const& dt,
	Obstacle *obstacle )
{
  list<MonObstacle*> allObstacles = obstacle->getObstacles();
  for (list<MonObstacle*>::iterator iv=allObstacles.begin();
  	iv!=allObstacles.end();iv++) (*iv)->setIndicator( 0. );
  obstacle->updateIndicator( temps, dt );
}




/* Ecriture des particules
--------------------------*/
void Paraview_PostProcessingWriter::writeParticulesPostProcessing_Paraview(
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
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0, i;
  for (particule=particules->begin();particule!=particules->end();particule++)
  {
    if ( (*particule)->getActivity() == COMPUTE &&
    	( (*particule)->getTag() != 2 || forceForAllTag ) )
    {
      nbpts += (*particule)->numberOfPoints_PARAVIEW();
      nbcells += (*particule)->numberOfCells_PARAVIEW();
    }
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;

  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point> ppp;
    list<Point>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
	( (*particule)->getTag() != 2 || forceForAllTag ) )
      {
        ppp = (*particule)->get_polygonsPts_PARAVIEW( PPTranslation );
        for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
          for (int comp=0;comp<3;++comp)
	    write_double_binary( (*ilpp)[comp] ) ;
      }
    flush_binary( f, "writeParticulesPostProcessing_Paraview/Points" );
  }
  else
    for (particule=particules->begin();particule!=particules->end();particule++)
      if ((*particule)->getActivity() == COMPUTE &&
	( (*particule)->getTag() != 2 || forceForAllTag ) )
        (*particule)->write_polygonsPts_PARAVIEW( f, PPTranslation );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity, offsets, cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber = 0, last_offset = 0;
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE &&
	( (*particule)->getTag() != 2 || forceForAllTag ) )
      (*particule)->write_polygonsStr_PARAVIEW(connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeParticulesPostProcessing_Paraview/connectivity" );
  }
  else
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeParticulesPostProcessing_Paraview/offsets" );
  }
  else
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeParticulesPostProcessing_Paraview/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;

  if( Grains_Exec::m_withFluidTemperature ||
      Grains_Exec::m_withSolidTemperature )
    f << "<CellData Scalars=\"NormU,NormOm,CoordNumb,Temperature\">" << endl;
  else
    f << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">" << endl;

  f << "<DataArray type=\"Float32\" Name=\"NormU\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE &&
       ( (*particule)->getTag() != 2 || forceForAllTag ) )
    {
      double normU = Norm( *(*particule)->getVitesseTranslation() );
      int nc = (*particule)->numberOfCells_PARAVIEW();
      if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( normU );
      else for (i=0;i<nc;++i) f << normU << " ";
    }
  if ( m_binary ) flush_binary( f,
  	"writeParticulesPostProcessing_Paraview/NormU" );
  f << endl;
  f << "</DataArray>" << endl;

  f << "<DataArray type=\"Float32\" Name=\"NormOm\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE &&
       ( (*particule)->getTag() != 2 || forceForAllTag ) )
    {
      double normOm = Norm( *(*particule)->getVitesseRotation() );
      int nc = (*particule)->numberOfCells_PARAVIEW();
      if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( normOm );
      else for (i=0;i<nc;++i) f << normOm << " ";
    }
  if( m_binary )
    flush_binary( f, "writeParticulesPostProcessing_Paraview/NormOm" );
  f << endl;
  f << "</DataArray>" << endl;

  f << "<DataArray type=\"Float32\" Name=\"CoordNumb\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE &&
       ( (*particule)->getTag() != 2 || forceForAllTag ) )
    {
      double coordNum = double((*particule)->getCoordinationNumber());
      int nc = (*particule)->numberOfCells_PARAVIEW();
      if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( coordNum );
      else for (i=0;i<nc;++i) f << coordNum << " ";
    }
  if ( m_binary ) flush_binary( f,
  	"writeParticulesPostProcessing_Paraview/CoordNumb" );
  f << endl;
  f << "</DataArray>" << endl;

  if( Grains_Exec::m_withFluidTemperature ||
      Grains_Exec::m_withSolidTemperature )
  {
    f << "<DataArray type=\"Float32\" Name=\"Temperature\" ";
    if( m_binary )
      f << "offset=\"" << OFFSET << "\" format=\"appended\">";
    else
      f << "format=\"ascii\">";
    f << endl;
    if( m_binary )
      start_output_binary( sizeof_Float32, int(cellstype.size()) );
    for (particule=particules->begin();particule!=particules->end();particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
      {
        double const* temperature = (*particule)->get_solidTemperature();
        int nc = (*particule)->numberOfCells_PARAVIEW();
        if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( *temperature );
        else for (i=0;i<nc;++i) f << *temperature << " ";
      }
    if( m_binary )
      flush_binary( f, "writeParticulesPostProcessing_Paraview/Temperature" );
    f << endl;
    f << "</DataArray>" << endl;
  }

  f << "</CellData>" << endl;

  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();
}




/* Ecriture du maillage de cellules du LinkedCell
-------------------------------------------------*/
void Paraview_PostProcessingWriter::writeLinkedCellPostProcessing_Paraview(
	LinkedCell const* LC, const string &partFilename )
{
  vector<Cellule*> const* allCells = LC->getAllCellules();
  vector<Cellule*>::const_iterator icell;
  Transform CelPosition;
  Convex* convexCellule = new Box( LC->getCelluleSize(X), LC->getCelluleSize(Y),
    	LC->getCelluleSize(Z) );
  Forme CelForme( convexCellule, CelPosition );
  Point const* cg = NULL;

  ofstream f( ( m_ParaviewFilename_root + "/" + partFilename ).c_str(),
  	ios::out );

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0;
  for (icell=allCells->begin();icell!=allCells->end();icell++)
  {
    nbpts+=8;
    nbcells+=1;
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point> ppp;
    list<Point>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (icell=allCells->begin();icell!=allCells->end();icell++)
    {
      cg = (*icell)->getCentre();
      CelForme.setOrigin( (*cg)[X], (*cg)[Y], (*cg)[Z] );
      ppp = CelForme.get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/Points" );
  }
  else
    for (icell=allCells->begin();icell!=allCells->end();icell++)
    {
      cg = (*icell)->getCentre();
      CelForme.setOrigin( (*cg)[X], (*cg)[Y], (*cg)[Z] );
      CelForme.write_polygonsPts_PARAVIEW( f );
    }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (icell=allCells->begin();icell!=allCells->end();icell++)
    CelForme.getConvex()->write_polygonsStr_PARAVIEW(
    	connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/connectivity" );
  }
  else
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/offsets" );
  }
  else
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();
}




/* Ecriture des fenetres d'insertion
------------------------------------*/
void Paraview_PostProcessingWriter::writeInsertionPostProcessing_Paraview(
   	vector<Fenetre> const& insert_windows,
  	const string &partFilename )
{
  vector<Fenetre>::const_iterator iv;
  list<Forme*> iwlist;
  list<Forme*>::iterator il = iwlist.begin();
  Transform gcwindow;
  Convex* ccw = NULL;
  Forme* ffw = NULL;
  Vecteur v_axis;
  Matrix mrot;
  Scalar Le = 0., Lh = 0., Lt = 0., angle = 0.,
  	mean_radius = 0., thickness = 0. ;
  Point panelCenter, center;
  size_t nbPanels = 32;
  Scalar polygonAngle = 2. * PI / Scalar(nbPanels);

  for (iv=insert_windows.begin();iv!=insert_windows.end();iv++)
  {
    switch( iv->ftype )
    {
      case FENETRE_BOX:
        gcwindow.setOrigin( 0.5 * ( iv->ptA + iv->ptB ) );
        ccw = new Box( 0.5 * ( iv->ptB - iv->ptA ) );
        ffw = new Forme( ccw, gcwindow );
        iwlist.push_back( ffw );
        break;

      case FENETRE_CYLINDER:
        switch ( iv->axisdir )
        {
          case X:
            v_axis[X] = iv->hauteur;
            mrot.setValue( cos( 0.5 * PI ), -sin( 0.5 * PI ), 0.,
                sin( 0.5 * PI ), cos( 0.5 * PI ), 0.,
                0., 0., 1. );
            break;

          case Y:
            v_axis[Y] = iv->hauteur;
            break;

          default:
            v_axis[Z] = iv->hauteur;
            mrot.setValue( 1., 0., 0.,
                0., cos( 0.5 * PI ), -sin( 0.5 * PI ),
                0., sin( 0.5 * PI ), cos( 0.5 * PI ) );
            break;
        }
        gcwindow.setOrigin( iv->ptA + 0.5 * v_axis );
        ccw = new Cylinder( iv->radius, iv->hauteur );
        ffw = new Forme( ccw, gcwindow );
        ffw->getTransform()->setBasis( mrot );
        iwlist.push_back( ffw );
        break;

      case FENETRE_ANNULUS:
        mean_radius = iv->radius_int + 0.5 * ( iv->radius - iv->radius_int ) ;
        thickness = iv->radius - iv->radius_int	;
        center = iv->ptA;
        Le = thickness;
        Lt = iv->radius * tan( polygonAngle );
        Lh = iv->hauteur;
        switch ( iv->axisdir )
        {
          case X:
            center[X] += 0.5 * iv->hauteur;
            for (size_t iNb=0; iNb!=nbPanels; iNb++)
            {
              angle = 2. * PI * Scalar(iNb) / Scalar(nbPanels);
              panelCenter[X] = center[X];
              panelCenter[Y] = center[Y] + mean_radius * cos(angle);
              panelCenter[Z] = center[Z] + mean_radius * sin(angle);
              gcwindow.setOrigin( panelCenter );

              mrot.setValue( 1., 0., 0.,
                  0., cos(angle), -sin(angle),
                  0., sin(angle), cos(angle) );
              gcwindow.setBasis( mrot );
              ccw = new Box( Lh, Le, Lt );
              ffw = new Forme( ccw, gcwindow );
              iwlist.push_back( ffw );
            }
            break;

          case Y:
            center[Y] += 0.5 * iv->hauteur;
            for (size_t iNb=0; iNb!=nbPanels; iNb++)
            {
              angle = 2. * PI * Scalar(iNb) / Scalar(nbPanels);
              panelCenter[X] = center[X] + mean_radius * sin(angle);
              panelCenter[Y] = center[Y];
              panelCenter[Z] = center[Z] + mean_radius * cos(angle);
              gcwindow.setOrigin( panelCenter );

              mrot.setValue( cos(angle), 0., sin(angle),
                  0., 1., 0.,
                  -sin(angle), 0., cos(angle) );
              gcwindow.setBasis( mrot );
              ccw = new Box( Lt, Lh, Le );
              ffw = new Forme( ccw, gcwindow );
              iwlist.push_back( ffw );
            }
            break;

          default:
            center[Z] += 0.5 * iv->hauteur;
            for (size_t iNb=0; iNb!=nbPanels; iNb++)
            {
              angle = 2. * PI * Scalar(iNb) / Scalar(nbPanels);
              panelCenter[X] = center[X] + mean_radius * cos(angle);
              panelCenter[Y] = center[Y] + mean_radius * sin(angle);
              panelCenter[Z] = center[Z];
              gcwindow.setOrigin( panelCenter );

              mrot.setValue( cos(angle), -sin(angle), 0.,
                  sin(angle), cos(angle), 0.,
                  0., 0., 1. );
              gcwindow.setBasis( mrot );
              ccw = new Box( Le, Lt, Lh );
              ffw = new Forme( ccw, gcwindow );
              iwlist.push_back( ffw );
            }
            break;
        }
        break;

      case FENETRE_LINE:
        gcwindow = Segment::computeTransform(
            0.5 * ( iv->ptB - iv->ptA ), 0.5 * ( iv->ptA + iv->ptB ) );
        ccw = new Segment( Norm( iv->ptB - iv->ptA ) );
        ffw = new Forme( ccw, gcwindow );
        iwlist.push_back( ffw );
        break;

      default:
        break;
    }
  }

  ofstream f( ( m_ParaviewFilename_root + "/" + partFilename + ".vtu" ).c_str(),
      ios::out );

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0;
  for (il=iwlist.begin();il!=iwlist.end();il++)
  {
    nbpts += (*il)->getConvex()->numberOfPoints_PARAVIEW();
    nbcells += (*il)->getConvex()->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point> ppp;
    list<Point>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (il=iwlist.begin();il!=iwlist.end();il++)
    {
      ppp = (*il)->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
          write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeInsertionPostProcessing_Paraview/Points" );
  }
  else
    for (il=iwlist.begin();il!=iwlist.end();il++)
      (*il)->write_polygonsPts_PARAVIEW( f );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (il=iwlist.begin();il!=iwlist.end();il++)
    (*il)->getConvex()->write_polygonsStr_PARAVIEW( connectivity,
        offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/connectivity" );
  }
  else
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/offsets" );
  }
  else
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();

  for (il=iwlist.begin();il!=iwlist.end();il++) delete *il;
}




/* Ecriture des vecteurs vitesse de translation & rotation des particules
-------------------------------------------------------------------------*/
void Paraview_PostProcessingWriter::
	writeVectorsMotionPostProcessing_Paraview(
	list<Particule*> const* particules, const string &partFilename )
{
  list<Particule*>::const_iterator particule;
  Point gc;
  Vecteur const* vectrans;
  Vecteur const* PPTranslation =
      Grains_Exec::m_translationParaviewPostProcessing ;

  ofstream f( ( m_ParaviewFilename_root + "/" + partFilename ).c_str(),
  	ios::out );

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0;
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE && (*particule)->getTag() != 2 )
      nbpts++;
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
        particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
           (*particule)->getTag() != 2 )
      {
        gc = *(*particule)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          write_double_binary( gc[comp] ) ;
      }
    flush_binary( f, "writeVectorsMotionPostProcessing_Paraview/Points" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         (*particule)->getTag() != 2 )
      {
        gc = *(*particule)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          f << gc[comp] << " " ;
      }
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 3 ) ;
      for (int ii=0;ii<3;ii++) write_int_binary( 0 );
      flush_binary( f,
      	"writeVectorsMotionPostProcessing_Paraview/connectivity" );
    }
    else
      f << "0 0 0";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1) ;
      write_int_binary( 3 );
      flush_binary( f, "writeVectorsMotionPostProcessing_Paraview/offsets" );
    }
    else
      f << "3";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 5 );
      flush_binary( f, "writeVectorsMotionPostProcessing_Paraview/types" );
    }
    else
      f << "5";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;

  f << "<PointData Vectors=\"U,Omega\">" << endl;
  f << "<DataArray Name=\"U\" NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE
		&& (*particule)->getTag() != 2 )
      {
        vectrans = (*particule)->getVitesseTranslation();
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vectrans)[comp] ) ;
      }
    flush_binary( f, "writeVectorsMotionPostProcessing_Paraview/U" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE
		&& (*particule)->getTag() != 2 )
      {
        vectrans = (*particule)->getVitesseTranslation();
        for (int comp=0;comp<3;++comp)
	  f << (*vectrans)[comp] << " " ;
      }
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "<DataArray Name=\"Omega\" NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE
		&& (*particule)->getTag() != 2 )
      {
        vectrans = (*particule)->getVitesseRotation();
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vectrans)[comp] ) ;
      }
    flush_binary( f, "writeVectorsMotionPostProcessing_Paraview/Omega" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE
		&& (*particule)->getTag() != 2 )
      {
        vectrans = (*particule)->getVitesseRotation();
        for (int comp=0;comp<3;++comp)
	  f << (*vectrans)[comp] << " " ;
      }
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "</PointData>" << endl;

  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();
}




/* Ecriture du reseau de forces de contact entre les particules
   J. F. Peters et al. Characterization of force chains in granular material.
   Phys. Rev. 2005
 * D. RAKOTONIRINA - Juin. 2017 - Creation
----------------------------------------------------------*/
void Paraview_PostProcessingWriter::writeForceChain_Paraview(
	list<Particule*> const* particules,
  	LinkedCell const* LC, const string &filename, Scalar const& temps,
    Scalar dt )
{
  bool binary = false;
  list<struct PointForcePostProcessing>* pallContacts =
  	LC->CalculerForcesPostProcessing( particules, dt );
  list<struct PointForcePostProcessing>::iterator contact;
  ofstream f( ( m_ParaviewFilename_root + "/" + filename ).c_str(),
  	ios::out );
  size_t nbContact = pallContacts->size();
  int nbContactsWithObstacles = 0;
  for (contact=pallContacts->begin();contact!=pallContacts->end();contact++){
    if (contact->contactWithObstacle) nbContactsWithObstacles++;
  }
  f << "<?xml version=\"1.0\"?>" << endl;
  f << "<VTKFile type=\"PolyData\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<PolyData>" << endl;
  f << "<Piece NumberOfPoints=\"" << 2*nbContact
    << "\" NumberOfVerts=\"" << 0
    << "\" NumberOfLines=\"" << nbContact
    << "\" NumberOfStrips=\"" << 0
    << "\" NumberOfPolys=\"" << 0
    << "\" NumberOfContactsWithObstacles=\"" << nbContactsWithObstacles <<
    "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\"";
  if ( binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << " format=\"ascii\">";
  f << endl;
  if ( binary ) start_output_binary( sizeof_Float32, 6*int(nbContact) );
  for (contact=pallContacts->begin();contact!=pallContacts->end();contact++)
  {
    // Point 1
    if ( contact->comp0->isObstacle() )
    {
      if ( binary )
	for ( int i=0; i<3; ++i )
	  write_double_binary( contact->geometricPointOfContact[i] );
      else
	for ( int i=0; i<3; ++i )
	  f << contact->geometricPointOfContact[i] << " " ;
    }
    else
    {
      if ( binary )
	for ( int i=0; i<3; ++i )
	  write_double_binary( (*contact->comp0->getPosition())[i] );
      else
	for ( int i=0; i<3; ++i )
	  f << (*contact->comp0->getPosition())[i] << " " ;
    }
    // Point 2
    if ( contact->comp1->isObstacle() )
    {
      if ( binary )
	for ( int i=0; i<3; ++i )
	  write_double_binary( contact->geometricPointOfContact[i] );
      else
	for ( int i=0; i<3; ++i )
	  f << contact->geometricPointOfContact[i] << " " ;
    }
    else
    {
      if ( binary )
	for ( int i=0; i<3; ++i )
	  write_double_binary( (*contact->comp1->getPosition())[i] );
      else
	for ( int i=0; i<3; ++i )
	  f << (*contact->comp1->getPosition())[i] << " " ;
    }
  }
  if ( binary ) flush_binary( f, "writeForceChain_Paraview/Points" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  f << "<PointData Scalars=\"Norm_F\">" << endl;
  f << "<DataArray type=\"Float32\" Name=\"Norm_F\" ";
  if ( binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;

  if ( binary ) start_output_binary( sizeof_Float32, 2*int(nbContact) );
  for (contact=pallContacts->begin();contact!=pallContacts->end();contact++)
  {
    if ( binary )
      for ( int i=0; i<2; ++i )
	write_double_binary( Norm( contact->contactForce ) );
    else
      for ( int i=0; i<2; ++i )
	f << Norm( contact->contactForce ) << " " ;
  }
  if ( binary ) flush_binary( f, "writeForceChain_Paraview/Norm_F" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</PointData>" << endl;

  f << "<Lines>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( binary )
  {
    start_output_binary( sizeof_Int32, 2*int(nbContact) );
    for ( int i=0; i<2*int(nbContact); i++ ) write_int_binary( i );
  }
  else for ( size_t i=0; i<2*nbContact; i++ ) f << i << " ";
  if ( binary ) flush_binary( f, "writeForceChain_Paraview/connectivity" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( binary )
  {
    start_output_binary( sizeof_Int32, int(nbContact) );
    for ( int i=0; i<int(nbContact); i++ ) write_int_binary( 2*i );
  }
  else for ( size_t i=1; i<=nbContact; ++i ) f << 2*i << " ";
  if ( binary ) flush_binary( f, "writeForceChain_Paraview/offsets" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Lines>" << endl;

  f << "</Piece>" << endl;
  f << "</PolyData>" << endl;
  if ( binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;

  f.close();
}




/* Ecriture des vecteurs force de contact entre composants
----------------------------------------------------------*/
void Paraview_PostProcessingWriter::writeVectorsForcePostProcessing_Paraview(
	list<Particule*> const* particules,
  	LinkedCell const* LC, const string &partFilename, Scalar const& temps,
    Scalar dt)
{
  list<struct PointForcePostProcessing>* pallContacts =
  	LC->CalculerForcesPostProcessing( particules, dt );
  list<struct PointForcePostProcessing>::iterator contact;
  Vecteur vec;
  Scalar total_dissip = 0.;

  ofstream f( ( m_ParaviewFilename_root + "/" + partFilename ).c_str(),
  	ios::out );

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = int(pallContacts->size()), nbcells = 0;
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (contact=pallContacts->begin();contact!=pallContacts->end();contact++)
    {
      for (int comp=0;comp<3;++comp)
	write_double_binary( contact->geometricPointOfContact[comp] ) ;
      total_dissip += contact->contactDissip;
    }
    flush_binary( f, "writeVectorsForcePostProcessing_Paraview/Points" );
  }
  else
  {
    for (contact=pallContacts->begin();contact!=pallContacts->end();contact++)
    {
      for (int comp=0;comp<3;++comp)
	 f << contact->geometricPointOfContact[comp] << " " ;
      total_dissip += contact->contactDissip;
    }
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 3 ) ;
      for (int ii=0;ii<3;ii++) write_int_binary( 0 );
      flush_binary( f,
      	"writeVectorsForcePostProcessing_Paraview/connectivity" );
    }
    else
      f << "0 0 0";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1) ;
      write_int_binary( 3 );
      flush_binary( f, "writeVectorsForcePostProcessing_Paraview/offsets" );
    }
    else
      f << "3";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 5 );
      flush_binary( f, "writeVectorsForcePostProcessing_Paraview/types" );
    }
    else
      f << "5";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "<PointData Vectors=\"Force\">" << endl;
  f << "<DataArray Name=\"Force\" NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (contact=pallContacts->begin();contact!=pallContacts->end();contact++)
      for (int comp=0;comp<3;++comp)
	write_double_binary( contact->contactForce[comp] ) ;
    flush_binary( f, "writeVectorsForcePostProcessing_Paraview/Force" );
  }
  else
  {
    for (contact=pallContacts->begin();contact!=pallContacts->end();contact++)
      for (int comp=0;comp<3;++comp)
        f << contact->contactForce[comp] << " " ;
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "</PointData>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();

  if ( Grains_Exec::m_ContactDissipation )
  {
    ofstream g( ( m_ParaviewFilename_root + "/ContEnerDiss_time.res" ).c_str(),
  	ios::app );
    g << temps << " " << total_dissip << endl;
    g.close();
  }
  pallContacts->clear();
  delete pallContacts;
}




/* Ecriture des particules de forme sph�rique sous forme d'un vecteur
---------------------------------------------------------------------*/
void Paraview_PostProcessingWriter:: writeSpheresPostProcessing_Paraview(
    list<Particule*> const* particules,
    const string &partFilename,
    bool const& forceForAllTag )
{
  list<Particule*>::const_iterator particule;
  Point gc;
  Vecteur vectrans;
  Vecteur const* PPTranslation =
  	Grains_Exec::m_translationParaviewPostProcessing ;

  ofstream f( ( m_ParaviewFilename_root + "/" + partFilename ).c_str(),
  	ios::out );

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0;
  for (particule=particules->begin();particule!=particules->end();particule++)
    if ( (*particule)->getActivity() == COMPUTE &&
       ( (*particule)->getTag() != 2 || forceForAllTag ) )
      nbpts++;
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
        particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
      {
        gc = *(*particule)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          write_double_binary( gc[comp] ) ;
      }
    flush_binary( f, "writeSpheresPostProcessing_Paraview/Points" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
        particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
           ( (*particule)->getTag() != 2 || forceForAllTag ) )
      {
        gc = *(*particule)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          f << gc[comp] << " " ;
      }
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 3 ) ;
      for (int ii=0;ii<3;ii++) write_int_binary( 0 );
      flush_binary( f, "writeSpheresPostProcessing_Paraview/connectivity" );
    }
    else
      f << "0 0 0";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 3 );
      flush_binary( f, "writeSpheresPostProcessing_Paraview/offsets" );
    }
    else
      f << "3";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 5 );
      flush_binary( f, "writeSpheresPostProcessing_Paraview/types" );
    }
    else
      f << "5";
  }
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "<PointData Vectors=\"Orientation\" ";
  if( Grains_Exec::m_withFluidTemperature ||
      Grains_Exec::m_withSolidTemperature )
    f << "Scalars=\"NormU,NormOm,CoordNumb,Temperature\">" << endl;
  else
    f << "Scalars=\"NormU,NormOm,CoordNumb\">" << endl;
  f << "<DataArray Name=\"Orientation\" "
    << "NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
      {
        vectrans = (*particule)->vecteurOrientation();
        for (int comp=0;comp<3;++comp)
          write_double_binary( vectrans[comp] ) ;
      }
    flush_binary( f, "writeSpheresPostProcessing_Paraview/Orientation" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
      {
        vectrans = (*particule)->vecteurOrientation();
        for (int comp=0;comp<3;++comp)
          f << vectrans[comp] << " " ;
      }
    f << endl;
  }
  f << "</DataArray>" << endl;

  f << "<DataArray Name=\"NormU\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
        particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
        write_double_binary( Norm( *(*particule)->getVitesseTranslation() ) );
    flush_binary( f, "writeSpheresPostProcessing_Paraview/NormU" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
        f << Norm( *(*particule)->getVitesseTranslation() ) << " " ;
    f << endl;
  }
  f << "</DataArray>" << endl;

  f << "<DataArray Name=\"NormOm\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
        write_double_binary( Norm( *(*particule)->getVitesseRotation() ) );
    flush_binary( f, "writeSpheresPostProcessing_Paraview/NormOm" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
        particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
        f << Norm( *(*particule)->getVitesseRotation() ) << " " ;
    f << endl;
  }
  f << "</DataArray>" << endl;

  f << "<DataArray Name=\"CoordNumb\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, nbpts ) ;
    for (particule=particules->begin();particule!=particules->end();
      	particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
        write_double_binary( double((*particule)->getCoordinationNumber()) );
    flush_binary( f, "writeSpheresPostProcessing_Paraview/CoordNumb" );
  }
  else
  {
    for (particule=particules->begin();particule!=particules->end();
        particule++)
      if ( (*particule)->getActivity() == COMPUTE &&
         ( (*particule)->getTag() != 2 || forceForAllTag ) )
        f << double((*particule)->getCoordinationNumber()) << " " ;
    f << endl;
  }
  f << "</DataArray>" << endl;

  if( Grains_Exec::m_withFluidTemperature ||
      Grains_Exec::m_withSolidTemperature )
  {
    f << "<DataArray Name=\"Temperature\" type=\"Float32\" ";
    if( m_binary )
      f << "offset=\"" << OFFSET << "\" format=\"appended\">";
    else
      f << "format=\"ascii\">";
    f << endl;
    if( m_binary )
    {
      start_output_binary( sizeof_Float32, nbpts ) ;
      for (particule=particules->begin();particule!=particules->end();
          particule++)
        if ( (*particule)->getActivity() == COMPUTE &&
           ( (*particule)->getTag() != 2 || forceForAllTag ) )
          write_double_binary( *((*particule)->get_solidTemperature()) );
      flush_binary( f, "writeSpheresPostProcessing_Paraview/Temperature" );
    }
    else
    {
      for (particule=particules->begin();particule!=particules->end();
          particule++)
        if ( (*particule)->getActivity() == COMPUTE &&
           ( (*particule)->getTag() != 2 || forceForAllTag ) )
          f << *((*particule)->get_solidTemperature()) << " " ;
      f << endl;
    }
    f << "</DataArray>" << endl;
  }

  f << "</PointData>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();
}




/* Ecriture du fichier pvtu
---------------------------*/
void Paraview_PostProcessingWriter::writePVTP_Paraview( const string &filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar )
{
  ofstream f( ( m_ParaviewFilename_root + "/" + filename + ".pvtp" ).c_str(),
  	ios::out );
  f << "<?xml version=\"1.0\"?>" << endl;
  f << "<VTKFile type=\"PPolyData\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<PPolyData GhostLevel=\"0\">" << endl;
  f << "<PPoints>" << endl;
  f << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << endl;
  f << "</PPoints>" << endl;
  f << "<PPointData Scalars=\"Norm_F\">" << endl;
  f << "<PDataArray type=\"Float32\" Name=\"Norm_F\"/>" << endl;
  f << "</PPointData>" << endl;
  f << "<PLines>" << endl;
  f << "<PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">"
  	<< endl;
  f << "</PDataArray>" << endl;
  f << "<PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">" << endl;
  f << "</PDataArray>" << endl;
  f << "</PLines>" << endl;
  for (int i=0;i<m_nprocs;++i)
  {
    // Does this processor have to write down outputs ?
    if( PostProcessingWriter::m_bPPWindow[i] )
    {
      f << "<Piece Source=\"" << filename << "_" << i << ".vtp\">" << endl;
      f << "</Piece>" << endl;
    }
  }
  f << "</PPolyData>" << endl;
  f << "</VTKFile>" << endl;
  f.close();
}




/* Ecriture du fichier pvtu
---------------------------*/
void Paraview_PostProcessingWriter::writePVTU_Paraview( const string &filename,
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
    // Does this processor have to write down outputs ?
    if( PostProcessingWriter::m_bPPWindow[i] )
    {
      f << "<Piece Source=\"" << filename << "_" << i << ".vtu\">" << endl;
      f << "</Piece>" << endl;
    }
  }
  f << "</PUnstructuredGrid>" << endl;
  f << "</VTKFile>" << endl;
  f.close();
}




/* Ecriture des composants pour une erreur de contact
-----------------------------------------------------*/
void Paraview_PostProcessingWriter::writeErreurComposantsPostProcessing(
	string const& filename,
  	list<Composant*> const& errcomposants )
{
  list<Composant*>::const_iterator compo;

  ostringstream ossRK;
  ossRK << m_rank;

  ofstream f( ( m_ParaviewFilename_root + "/" + filename + "_" + ossRK.str()
  	+ ".vtu" ).c_str(), ios::out );
  f << "<?xml version=\"1.0\"?>" << endl;
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0;
  for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
  {
//    nbpts += (*compo)->getForme()->getConvex()->numberOfPoints_PARAVIEW();
//    nbcells += (*compo)->getForme()->getConvex()->numberOfCells_PARAVIEW();
    nbpts += (*compo)->numberOfPoints_PARAVIEW();
    nbcells += (*compo)->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point> ppp;
    list<Point>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
    {
      //ppp = (*compo)->getForme()->get_polygonsPts_PARAVIEW();
      ppp = (*compo)->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeErreurComposantsPostProcessing/Points" );
  }
  else
    for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
      (*compo)->write_polygonsPts_PARAVIEW(f);
      //(*compo)->getForme()->write_polygonsPts_PARAVIEW(f);
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
    (*compo)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );
  //  (*compo)->getForme()->getConvex()->write_polygonsStr_PARAVIEW(
  //  	connectivity, offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeErreurComposantsPostProcessing/connectivity" );
  }
  else
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeErreurComposantsPostProcessing/offsets" );
  }
  else
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeErreurComposantsPostProcessing/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;
  }
  f << "</VTKFile>" << endl;
  f.close();
}




/* Efface les fichiers resultats
--------------------------------*/
void Paraview_PostProcessingWriter::clearResultFiles() const
{
  if ( m_rank == 0 )
  {
    string cmd = "bash " + Grains_Exec::m_GRAINS_HOME
     	+ "/ExecScripts/Paraview_clear.exec " + m_ParaviewFilename_root +
	" " + m_ParaviewFilename;
    system( cmd.c_str() );
  }
}




/* Numero de cycle initial
--------------------------*/
void Paraview_PostProcessingWriter::setInitialCycleNumber( const int& cycle0 )
{
  m_ParaviewCycleNumber = cycle0;
  m_initialCycleNumber_forced = true ;
}




/* Ecriture dans un fichier d'un long ostringstream (routine
   developpee suite a des prob de taille de buffer sur ENER110)
---------------------------------------------------------------*/
void Paraview_PostProcessingWriter::writeBigOSS( ofstream& fileOUT,
	ostringstream const& oss )
{
  string str = oss.str();
  for ( string::iterator it=str.begin(); it!=str.end(); ++it)
    fileOUT << *it;
}




//-----------------------------------------------------------------------
// Methodes copiees de Pelicans-2.2.4: classe Paraview_PostProcessingWriter
//-----------------------------------------------------------------------
void Paraview_PostProcessingWriter:: start_output_binary( int size, int number )
{
  int current_output_size = size*number ;
//   unsigned long ncomp = current_output_size + (current_output_size+999)/1000
//   	+ 12 + sizeof_Int32 ;
  int ncomp = current_output_size + (current_output_size+999)/1000
  	+ 12 + sizeof_Int32 ;
  check_allocated_binary( ncomp ) ;
  CURRENT_LENGTH = store_int_binary(current_output_size) ;
}




void Paraview_PostProcessingWriter:: write_double_binary( double val )
{
  *((float*)&(BUFFER[OFFSET])) = (float)val ;
  OFFSET += sizeof_Float32  ;
}




void Paraview_PostProcessingWriter:: write_int_binary( int val )
{
  store_int_binary(val) ;
}




int Paraview_PostProcessingWriter:: store_int_binary( int val )
{
  int result = OFFSET ;
  *((int*)&(BUFFER[OFFSET])) = val ;
  OFFSET += sizeof_Int32  ;
  return result ;
}




void Paraview_PostProcessingWriter:: check_allocated_binary( int size )
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




void Paraview_PostProcessingWriter:: flush_binary( std::ofstream& file,
	string const& calling )
{
  compress_segment_binary(CURRENT_LENGTH,calling) ;
  file << endl ;
}




void Paraview_PostProcessingWriter:: compress_segment_binary( int seg,
	string const& calling )
{
   static int BlockSize = 32768 ;
   int size = (int)(*((int*)&BUFFER[seg])) ;

   int numFullBlocks = size / BlockSize;
   int lastBlockSize = size % BlockSize;
   int numBlocks = numFullBlocks + (lastBlockSize?1:0);

   int headerLength = numBlocks+3;

   int * CompressionHeader = new int[headerLength];
   CompressionHeader[0] = numBlocks;
   CompressionHeader[1] = BlockSize;
   CompressionHeader[2] = lastBlockSize;

   unsigned long encoded_buff_size = max(BlockSize,size)  ;
   unsigned char* encoded_buff = new unsigned char [ encoded_buff_size ] ;
   int encoded_offset = 0 ;
   for( int block=0 ; block<numBlocks ; block++ )
   {
      int buffer_start = seg + sizeof_Int32 + block*BlockSize ;
      int length = ( block+1<numBlocks || !lastBlockSize ?
      	BlockSize : lastBlockSize ) ;
      unsigned char* to_encode = (unsigned char *)(&BUFFER[buffer_start]) ;
      unsigned char* encoded = &encoded_buff[encoded_offset] ;
      unsigned long ncomp = encoded_buff_size - encoded_offset ;

      if(compress2((Bytef*)encoded,
                   &ncomp,
                   (const Bytef*)to_encode,
                   length,
                   Z_DEFAULT_COMPRESSION) != Z_OK)
      {
         cout << "Zlib error while compressing data." << endl;
	 cout << "from " << calling << endl;
	 cout << "Details : block = " << block << "  numBlocks = " <<
	 	numBlocks << "  length = " << length << "  ncomp = "
		<< ncomp << endl;
	 exit(0);
      }
//       CompressionHeader[3+block] = ncomp ;
//       encoded_offset += ncomp ;
      CompressionHeader[3+block] = int(ncomp) ;
      encoded_offset += int(ncomp) ;
   }

   OFFSET = seg ;
   check_allocated_binary( headerLength*sizeof_Int32 + encoded_offset ) ;

   for(int i=0 ; i<headerLength ; i++ )
      store_int_binary(CompressionHeader[i]) ;

   for(int i=0 ; i<encoded_offset ; i++ )
      BUFFER[OFFSET++] = encoded_buff[i] ;

   if( OFFSET%4 != 0 )
      OFFSET = 4*( OFFSET/4 +1 ) ; // Re-allignement

   delete [] CompressionHeader ;
   delete [] encoded_buff ;
}
