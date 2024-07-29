#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "Particle.hh"
#include "Component.hh"
#include "Obstacle.hh"
#include "Box.hh"
#include "Cylinder.hh"
#include "Segment.hh"
#include "Cell.hh"
#include "Vector3.hh"
#include <zlib.h>
using namespace solid;

static int sizeof_Float32 = 4 ;
static int sizeof_Int32 = 4 ;


// ----------------------------------------------------------------------------
// Writes particles data
void ParaviewPostProcessingWriter::writeParticles_Paraview(
	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag, bool const& processwrites )
{
  list<Particle*>::const_iterator particle;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;

  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename ).c_str(), 
  	ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0, i;
  for (particle=particles->begin();particle!=particles->end();particle++)
  {
    if ( (*particle)->getActivity() == COMPUTE && 
    	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nbpts += (*particle)->numberOfPoints_PARAVIEW();
      nbcells += (*particle)->numberOfCells_PARAVIEW();
    }
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;

  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    list<Point3> ppp;
    list<Point3>::iterator ilpp;    
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
      {
        ppp = (*particle)->get_polygonsPts_PARAVIEW( PPTranslation );
        for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
          for (int comp=0;comp<3;++comp)
	    write_double_binary( (*ilpp)[comp] ) ;
      }
    flush_binary( f, "writeParticles_Paraview/Points" );      
  }
  else
    for (particle=particles->begin();particle!=particles->end();particle++)
      if ((*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
        (*particle)->write_polygonsPts_PARAVIEW( f, PPTranslation );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity, offsets, cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber = 0, last_offset = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE &&
	( (*particle)->getTag() != 2 || forceForAllTag ) )    
      (*particle)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeParticles_Paraview/connectivity" );
  }
  else 
  { 
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";	
    f << endl; 
  }     
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeParticles_Paraview/offsets" );
  }
  else
  {  
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";	
    f << endl;
  } 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeParticles_Paraview/types" );
  }
  else 
  { 
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";	
    f << endl;
  }
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;

  f << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">" << endl;

  f << "<DataArray type=\"Float32\" Name=\"NormU\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">"; 
  else f << "format=\"ascii\">" << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      double normU = Norm( *(*particle)->getTranslationalVelocity() );
      int nc = (*particle)->numberOfCells_PARAVIEW();
      if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( normU );
      else for (i=0;i<nc;++i) f << normU << " ";
    }
  if ( m_binary ) flush_binary( f, 
  	"writeParticles_Paraview/NormU" ); 
  else f << endl;
  f << "</DataArray>" << endl;      

  f << "<DataArray type=\"Float32\" Name=\"NormOm\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">"; 
  else f << "format=\"ascii\">" << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      double normOm = Norm( *(*particle)->getAngularVelocity() );
      int nc = (*particle)->numberOfCells_PARAVIEW();
      if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( normOm );
      else for (i=0;i<nc;++i) f << normOm << " ";
    }
  if( m_binary )
    flush_binary( f, "writeParticles_Paraview/NormOm" ); 
  else f << endl;
  f << "</DataArray>" << endl; 

  f << "<DataArray type=\"Float32\" Name=\"CoordNumb\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">"; 
  else f << "format=\"ascii\">" << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      double coordNum = double((*particle)->getCoordinationNumber());
      int nc = (*particle)->numberOfCells_PARAVIEW();
      if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( coordNum );
      else for (i=0;i<nc;++i) f << coordNum << " ";
    }
  if ( m_binary ) flush_binary( f, 
  	"writeParticles_Paraview/CoordNumb" ); 
  else f << endl;
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




// ----------------------------------------------------------------------------
// Writes spherical particles data in a vector form containing the
// center of mass coordinates of each particle
void ParaviewPostProcessingWriter:: writeSpheres_Paraview(
    list<Particle*> const* particles,
    string const& partFilename,
    bool const& forceForAllTag, bool const& processwrites )
{
  list<Particle*>::const_iterator particle;
  Point3 gc; 
  Vector3 vectrans;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;

  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename ).c_str(), 
  	ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
      nbpts++;
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
        particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
      {
        gc = *(*particle)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          write_double_binary( gc[comp] ) ;
      }
    flush_binary( f, "writeSpheres_Paraview/Points" );
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
        particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
           ( (*particle)->getTag() != 2 || forceForAllTag ) )
      {
        gc = *(*particle)->getPosition();
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
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 3 ) ;
      for (int ii=0;ii<3;ii++) write_int_binary( 0 );
      flush_binary( f, "writeSpheres_Paraview/connectivity" );
    }
    else  
      f << "0 0 0";
  }
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 3 );
      flush_binary( f, "writeSpheres_Paraview/offsets" );
    }
    else  
      f << "3";	
  }
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 5 );
      flush_binary( f, "writeSpheres_Paraview/types" );
    }
    else  
      f << "5";
  }
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "<PointData Vectors=\"Orientation\" ";
  f << "Scalars=\"NormU,NormOm,CoordNumb\">" << endl;
  f << "<DataArray Name=\"Orientation\" "
    << "NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
      {
        vectrans = (*particle)->computeOrientationVector();
        for (int comp=0;comp<3;++comp)
          write_double_binary( vectrans[comp] ) ;
      }
    flush_binary( f, "writeSpheres_Paraview/Orientation" );      
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
      {
        vectrans = (*particle)->computeOrientationVector();
        for (int comp=0;comp<3;++comp)
          f << vectrans[comp] << " " ;
        f << endl;	
      }
  }  
  f << "</DataArray>" << endl;

  f << "<DataArray Name=\"NormU\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
        particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
        write_double_binary( Norm( *(*particle)->getTranslationalVelocity() ) );
    flush_binary( f, "writeSpheres_Paraview/NormU" );      
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
        f << Norm( *(*particle)->getTranslationalVelocity() ) << " " ;
    f << endl;
  }
  f << "</DataArray>" << endl;

  f << "<DataArray Name=\"NormOm\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
        write_double_binary( Norm( *(*particle)->getAngularVelocity() ) );
    flush_binary( f, "writeSpheres_Paraview/NormOm" );
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
        particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
        f << Norm( *(*particle)->getAngularVelocity() ) << " " ;
    f << endl;
  }
  f << "</DataArray>" << endl;

  f << "<DataArray Name=\"CoordNumb\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
        write_double_binary( double((*particle)->getCoordinationNumber()) );
    flush_binary( f, "writeSpheres_Paraview/CoordNumb" );
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
        particle++)
      if ( (*particle)->getActivity() == COMPUTE && 
         ( (*particle)->getTag() != 2 || forceForAllTag ) )
        f << double((*particle)->getCoordinationNumber()) << " " ;
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




// ----------------------------------------------------------------------------
// Writes particle translational and angular velocity vectors
void ParaviewPostProcessingWriter::
	writeParticleVelocityVectors_Paraview(
	list<Particle*> const* particles, string const& partFilename,
	bool const& processwrites )
{
  list<Particle*>::const_iterator particle;
  Point3 gc; 
  Vector3 const* vec;
  Vector3 const* PPTranslation = 
      GrainsExec::m_translationParaviewPostProcessing ;

  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename ).c_str(), 
  	ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && (*particle)->getTag() != 2 )
      nbpts++;
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
        particle++)
      if ( (*particle)->getActivity() == COMPUTE &&
           (*particle)->getTag() != 2 )
      {
        gc = *(*particle)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          write_double_binary( gc[comp] ) ;
      }
    flush_binary( f, 
    	"writeParticleVelocityVectors_Paraview/Points" );
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE &&
         (*particle)->getTag() != 2 )
      {
        gc = *(*particle)->getPosition();
        if ( PPTranslation ) gc += *PPTranslation ;
        for (int comp=0;comp<3;++comp)
          f << gc[comp] << " " ;
        f << endl;	
      }
  }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 3 ) ;
      for (int ii=0;ii<3;ii++) write_int_binary( 0 );
      flush_binary( f, 
      	"writeParticleVelocityVectors_Paraview/connectivity" );
    }
    else  
      f << "0 0 0";
  }   
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1) ;
      write_int_binary( 3 );
      flush_binary( f, 
      	"writeParticleVelocityVectors_Paraview/offsets" );
    }
    else  
      f << "3";	
  }
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 5 );
      flush_binary( f, 
      	"writeParticleVelocityVectors_Paraview/types" );
    }
    else  
      f << "5";
  }
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  
  f << "<PointData Vectors=\"U,Omega\">" << endl;
  f << "<DataArray Name=\"U\" NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE 
		&& (*particle)->getTag() != 2 )
      {
        vec = (*particle)->getTranslationalVelocity();
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vec)[comp] ) ;
      }
    flush_binary( f, "writeParticleVelocityVectors_Paraview/U" );
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE 
		&& (*particle)->getTag() != 2 )
      {
        vec = (*particle)->getTranslationalVelocity();
        for (int comp=0;comp<3;++comp)
	  f << (*vec)[comp] << " " ;      
      }
    f << endl;
  }  
  f << "</DataArray>" << endl;
  f << "<DataArray Name=\"Omega\" NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE 
		&& (*particle)->getTag() != 2 )
      {
        vec = (*particle)->getAngularVelocity();
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*vec)[comp] ) ;
      }
    flush_binary( f, 
    	"writeParticleVelocityVectors_Paraview/Omega" );      
  }
  else
  {
    for (particle=particles->begin();particle!=particles->end();
      	particle++)
      if ( (*particle)->getActivity() == COMPUTE 
		&& (*particle)->getTag() != 2 )
      {
        vec = (*particle)->getAngularVelocity();
        for (int comp=0;comp<3;++comp)
	  f << (*vec)[comp] << " " ;      
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




// ----------------------------------------------------------------------------
// Writes contact force vectors
void ParaviewPostProcessingWriter::
	writeContactForceVectors_Paraview(
	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites )
{
  vector<struct PointForcePostProcessing> const* pallContacts = 
  	LC->getPPForces();
  size_t i = 0, comp = 0, nPPF = LC->getNbPPForces();
  Vector3 const* PPTranslation = 
      GrainsExec::m_translationParaviewPostProcessing ;
  Point3 pt;
          
  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename ).c_str(), 
  	ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  size_t nbpts = nPPF, nbcells = 0;
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">" << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3 * int(nbpts) ) ;
    for (i=0;i<nPPF;++i)
    {
      pt = (*pallContacts)[i].geometricPointOfContact;
      if ( PPTranslation ) pt += *PPTranslation ;      
      for (comp=0;comp<3;++comp)
	write_double_binary( pt[comp] ) ;
    }
    flush_binary( f, "writeContactForceVectors_Paraview/Points" );
  }
  else
  {
    for (i=0;i<nPPF;++i)
    {
      pt = (*pallContacts)[i].geometricPointOfContact;
      if ( PPTranslation ) pt += *PPTranslation ;       
      for (comp=0;comp<3;++comp)
	 f << pt[comp] << " " ; 
      f << endl;	 
    }	 
  }
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 3 ) ;
      for (size_t ii=0;ii<3;ii++) write_int_binary( 0 );
      flush_binary( f, 
      	"writeContactForceVectors_Paraview/connectivity" );
    }
    else  
      f << "0 0 0";
  }      
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1) ;
      write_int_binary( 3 );
      flush_binary( f, 
      	"writeContactForceVectors_Paraview/offsets" );
    }
    else  
      f << "3";	
  }
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  if ( nbpts )
  {
    if ( m_binary )
    {
      start_output_binary( sizeof_Int32, 1 ) ;
      write_int_binary( 5 );
      flush_binary( f, 
      	"writeContactForceVectors_Paraview/types" );
    }
    else  
      f << "5";
  }
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "<PointData Vectors=\"Force\">" << endl;
  f << "<DataArray Name=\"Force\" NumberOfComponents=\"3\" type=\"Float32\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;  
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3 * int(nbpts) ) ;   
    for (i=0;i<nPPF;++i)
      for (comp=0;comp<3;++comp)
	write_double_binary( (*pallContacts)[i].contactForce[comp] ) ;
    flush_binary( f, "writeContactForceVectors_Paraview/Force" );
  }
  else
  {
    for (i=0;i<nPPF;++i)
    {
      for (comp=0;comp<3;++comp)
        f << (*pallContacts)[i].contactForce[comp] << " " ;
      f << endl;
    }
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




// ----------------------------------------------------------------------------
// Writes force chain network data
void ParaviewPostProcessingWriter::writeContactForceChains_Paraview(
	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& filename, double const& time )
{
  vector<struct PointForcePostProcessing> const* pallContacts = 
  	LC->getPPForces();
  size_t i, k, nPPF = LC->getNbPPForces();
  Vector3 const* PPTranslation = 
      GrainsExec::m_translationParaviewPostProcessing ;
  Point3 pt0, pt1;
  int j;
  double norm = 0.;

  ofstream f( ( m_ParaviewFilename_dir + "/" + filename ).c_str(), 
  	ios::out );
  f << "<VTKFile type=\"PolyData\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<PolyData>" << endl;
  f << "<Piece NumberOfPoints=\"" << 2*nPPF 
    << "\" NumberOfVerts=\"" << 0
    << "\" NumberOfLines=\"" << nPPF
    << "\" NumberOfStrips=\"" << 0
    << "\" NumberOfPolys=\"" << 0 << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\"";
  if ( m_binary ) f << " offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << " format=\"ascii\">";
  f << endl;  
  if ( m_binary ) start_output_binary( sizeof_Float32, 6*int(nPPF) );
  for (i=0;i<nPPF;++i)
  {
    // Point 0
    pt0 = (*pallContacts)[i].PPptComp0;
    if ( PPTranslation ) pt0 += *PPTranslation ; 
    if ( m_binary )
      for (k=0;k<3;++k) write_double_binary( pt0[k] );
    else
      for (k=0;k<3;++k) f << pt0[k] << " " ;    

    // Point 1
    pt1 = (*pallContacts)[i].PPptComp1;
    if ( PPTranslation ) pt1 += *PPTranslation ; 
    if ( m_binary )
      for (k=0;k<3;++k) write_double_binary( pt1[k] );
    else
      for (k=0;k<3;++k) f << pt1[k] << " " ;
  }
  if ( m_binary ) flush_binary( f, "writeContactForceChains_Paraview/Points" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  f << "<PointData Scalars=\"ForceMag\">" << endl;
  f << "<DataArray type=\"Float32\" Name=\"ForceMag\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;

  if ( m_binary ) start_output_binary( sizeof_Float32, 2*int(nPPF) );
  for (i=0;i<nPPF;++i)
  {
    norm = Norm( (*pallContacts)[i].contactForce );
    if ( m_binary )
      for (k=0;k<2;++k)
	write_double_binary( norm );
    else
      for (k=0;k<2;++k)
	f << norm << " " ;
  }
  if ( m_binary ) flush_binary( f, 
  	"writeContactForceChains_Paraview/ForceMag" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</PointData>" << endl;

  f << "<Lines>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, 2*int(nPPF) );
    for (j=0;j<2*int(nPPF);j++) write_int_binary( j );
  }
  else 
    for (j=0;j<2*int(nPPF);j++) f << j << " ";
  if ( m_binary ) flush_binary( f, 
  	"writeContactForceChains_Paraview/connectivity" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(nPPF) );
    for (j=0;j<int(nPPF);j++) write_int_binary( 2 * j );
  }
  else for (j=0;j<int(nPPF);j++) f << 2*j << " ";
  if ( m_binary ) flush_binary( f, "writeContactForceChains_Paraview/offsets" );
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Lines>" << endl;

  f << "</Piece>" << endl;
  f << "</PolyData>" << endl;
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
