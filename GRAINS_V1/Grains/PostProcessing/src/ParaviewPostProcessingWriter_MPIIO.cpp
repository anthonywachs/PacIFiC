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
// Writes particles data in a single MPI file in text mode with 
// MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticlesPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  MPI_Datatype num_as_string;
  int nbpts = 0, nbcells = 0, nparts = 0, i, rank, j, nc, coordNum;
  list<Particle*>::const_iterator particle;
  char fmt[7] = "%8.6f ";
  char endfmt[7] = "%8.6f\n";
  const int charspernum = 9;
  size_t counter = 0;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  double nu, nom;     
  
  // Create a MPI datatype to write doublea as stringa with a given format
  MPI_Type_contiguous( charspernum, MPI_CHAR, &num_as_string ); 
  MPI_Type_commit( &num_as_string );

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points and cells
  for (particle=particles->begin();particle!=particles->end();particle++)
  {
    if ( (*particle)->getActivity() == COMPUTE && 
    	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nbpts += (*particle)->numberOfPoints_PARAVIEW();
      nbcells += (*particle)->numberOfCells_PARAVIEW();
      ++nparts;
    }
  }
  
  int total_nbpts = wrapper->sum_INT( nbpts ); 
  int total_nbcells = wrapper->sum_INT( nbcells );
  int* nbpts_per_proc = wrapper->AllGather_INT( nbpts ); 
  int* nparts_per_proc = wrapper->AllGather_INT( nparts );
  int* nbcells_per_proc = wrapper->AllGather_INT( nbcells );       


  // File header
  ostringstream oss;
  oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">\n";
  oss << "<UnstructuredGrid>\n";
  oss << "<Piece NumberOfPoints=\"" << total_nbpts << "\""
    	<< " NumberOfCells=\"" << total_nbcells << "\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "format=\"ascii\">";
  oss << "\n";        
  int header = int(oss.str().size());
  if ( m_rank == 0 )
    MPI_File_write( file, oss.str().c_str(), header, MPI_CHAR, &status );

  // Write point coordinates to the MPI file
  char* pts_coord = new char[ 3 * nbpts * charspernum + 1 ];
  list<Point3> ppp;
  list<Point3>::iterator ilpp;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      ppp = (*particle)->get_polygonsPts_PARAVIEW( PPTranslation );
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++,counter++)
      {
        sprintf( &pts_coord[3*counter*charspernum], fmt, (*ilpp)[X] );
	sprintf( &pts_coord[(3*counter+1)*charspernum], fmt, (*ilpp)[Y] );
	sprintf( &pts_coord[(3*counter+2)*charspernum], endfmt, (*ilpp)[Z] );
      }
    }    

  vector<int> mpifile_offsets( m_nprocs, 0 );
  mpifile_offsets[0] = header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], pts_coord, 3 * nbpts, 
    	num_as_string, &status ); 

  delete [] pts_coord;  


  // Header for connectivity
  ostringstream oss2;
  oss2 << "</DataArray>\n";
  oss2 << "</Points>\n";
  oss2 << "<Cells>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  oss2 << "format=\"ascii\">\n";
  header = int(oss2.str().size());
  int mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss2.str().c_str(), header, 
    	MPI_CHAR, &status );
    
  // Connectivity
  list<int> connectivity, offsets, cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber = 0, last_offset = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE &&
	( (*particle)->getTag() != 2 || forceForAllTag ) )    
      (*particle)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );  

  // Compute the point number shifts 
  vector<int> shift( m_nprocs, 0 );
  for (rank=1;rank<m_nprocs;rank++)
    for (j=0;j<rank;j++)
      shift[rank] += nbpts_per_proc[j];
  
  // Write connectivity to the MPI file with the point number shifts
  ostringstream* oss_out = new ostringstream;
  for (ii=connectivity.begin();ii!=connectivity.end();ii++)
    *oss_out << *ii + shift[m_rank] << " ";
  if ( m_rank == m_nprocs - 1 ) *oss_out << "\n";  
  int out_length = int(oss_out->str().size());
  int* out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  delete out_length_per_proc;


  // Header for offset
  ostringstream oss3;
  oss3 << "</DataArray>\n";
  oss3 << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  oss3 << "format=\"ascii\">\n";
  header = int(oss3.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ out_length_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss3.str().c_str(), header, 
    	MPI_CHAR, &status );

  // Compute the offset shifts 
  int max_offset = offsets.empty() ? 0 : offsets.back();
  int* max_offset_per_proc = wrapper->AllGather_INT( max_offset );
  for (rank=1;rank<m_nprocs;rank++)
  {
    shift[rank] = 0;
    for (j=0;j<rank;j++)
      shift[rank] += max_offset_per_proc[j];
  }  

  // Write offsets to the MPI file with the offset shifts
  oss_out = new ostringstream;
  for (ii=offsets.begin();ii!=offsets.end();ii++)
    *oss_out << *ii + shift[m_rank] << " ";
  if ( m_rank == m_nprocs - 1 ) *oss_out << "\n";  
  out_length = int(oss_out->str().size());
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  delete max_offset_per_proc;
  delete out_length_per_proc;  


  // Header for cell types
  ostringstream oss4;
  oss4 << "</DataArray>\n";
  oss4 << "<DataArray type=\"Int32\" Name=\"types\" ";
  oss4 << "format=\"ascii\">\n";
  header = int(oss4.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ out_length_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss4.str().c_str(), header, 
    	MPI_CHAR, &status );

  // Write cell types to the MPI file
  oss_out = new ostringstream;
  for (ii=cellstype.begin();ii!=cellstype.end();ii++)
    *oss_out << *ii << " ";
  if ( m_rank == m_nprocs - 1 ) *oss_out << "\n";  
  out_length = int(oss_out->str().size());
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  delete out_length_per_proc;


  // Field values on cells
  ostringstream oss5;
  oss5 << "</DataArray>\n";
  oss5 << "</Cells>\n";
  oss5 << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">\n";

  // Norm of translational velocity
  oss5 << "<DataArray type=\"Float32\" Name=\"NormU\" ";        
  oss5 << "format=\"ascii\">\n";    
  header = int(oss5.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ out_length_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss5.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  char* normU = new char[ nbcells * charspernum + 1 ];
  counter = 0;
  for (particle=particles->begin();
    	particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nc = (*particle)->numberOfCells_PARAVIEW();
      nu = Norm( *(*particle)->getTranslationalVelocity() );
      for (j=0;j<nc;++j)
      {
        sprintf( &normU[counter*charspernum], fmt, nu );
        ++counter;	
      }
    }

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbcells_per_proc[i-1] * charspernum * sizeof(char) ;
 
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], normU, nbcells, 
    	num_as_string, &status ); 

  delete [] normU; 
  
  
  // Norm of angular velocity
  ostringstream oss6;
  oss6 << "\n</DataArray>\n";  
  oss6 << "<DataArray type=\"Float32\" Name=\"NormOm\" ";        
  oss6 << "format=\"ascii\">\n";    
  header = int(oss6.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbcells_per_proc[m_nprocs-1] * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss6.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  char* normOm = new char[ nbcells * charspernum + 1 ];
  counter = 0;
  for (particle=particles->begin();
    	particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nc = (*particle)->numberOfCells_PARAVIEW();
      nom = Norm( *(*particle)->getAngularVelocity() );
      for (j=0;j<nc;++j)
      {
        sprintf( &normOm[counter*charspernum], fmt, nom );
        ++counter;	
      }
    }

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbcells_per_proc[i-1] * charspernum * sizeof(char) ;
 
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], normOm, nbcells, 
    	num_as_string, &status ); 

  delete [] normOm;
  
  
  // Coordination number
  ostringstream oss7;
  oss7 << "\n</DataArray>\n";  
  oss7 << "<DataArray type=\"Float32\" Name=\"CoordNumb\" ";        
  oss7 << "format=\"ascii\">\n";    
  header = int(oss7.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbcells_per_proc[m_nprocs-1] * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss7.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  oss_out = new ostringstream;
  for (particle=particles->begin();
    	particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nc = (*particle)->numberOfCells_PARAVIEW();
      coordNum = (*particle)->getCoordinationNumber();
      for (j=0;j<nc;++j)
        *oss_out << coordNum << " ";
    }
  out_length = int(oss_out->str().size());
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  delete out_length_per_proc;
   
   
  // Closing text
  ostringstream oss8;
  oss8 << "\n</DataArray>\n";
  oss8 << "</CellData>\n";
  oss8 << "</Piece>\n";
  oss8 << "</UnstructuredGrid>\n";
  oss8 << "</VTKFile>\n";
  header = int(oss8.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ out_length_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss8.str().c_str(), header, 
    	MPI_CHAR, &status );   


  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  MPI_Type_free( &num_as_string );
  delete [] nbpts_per_proc;
  delete [] nparts_per_proc; 
  delete [] nbcells_per_proc;    
}




// ----------------------------------------------------------------------------
// Writes particles data in a single MPI file in binary mode with 
// MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticlesPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  int nbpts = 0, nbcells = 0, nparts = 0, i, rank, nc;
  list<Particle*>::const_iterator particle;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  double normU, normOm, coordNum;
  int point_binary_offset = 0, connectivity_binary_offset, 
  	offsets_binary_offset, cellstype_binary_offset, normU_binary_offset,
	normOm_binary_offset, coord_binary_offset, total_offset;     

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points and cells
  for (particle=particles->begin();particle!=particles->end();particle++)
  {
    if ( (*particle)->getActivity() == COMPUTE && 
    	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nbpts += (*particle)->numberOfPoints_PARAVIEW();
      nbcells += (*particle)->numberOfCells_PARAVIEW();
      ++nparts;
    }
  }  


  // Write point coordinates to the binary buffer
  list<Point3> ppp;
  list<Point3>::iterator ilpp;    
  start_output_binary( sizeof_Float32, 3 * nbpts ) ;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      ppp = (*particle)->get_polygonsPts_PARAVIEW( PPTranslation );
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*ilpp)[comp] ) ;
    }
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticlesPostProcessing_Paraview_MPIIO_binary/Points" );


  // Write Connectivity to the binary buffer
  list<int> connectivity, offsets, cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber = 0, last_offset = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE &&
	( (*particle)->getTag() != 2 || forceForAllTag ) )    
      (*particle)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );  
	
  connectivity_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
  for (ii=connectivity.begin();ii!=connectivity.end();ii++)
    write_int_binary( *ii );	
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticlesPostProcessing_Paraview_MPIIO_binary/connectivity" );

  // Write offsets to the binary buffer  
  offsets_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, int(offsets.size()) ) ;  
  for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );  
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticlesPostProcessing_Paraview_MPIIO_binary/offsets" ); 


  // Write cell types to the binary buffer  
  cellstype_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;  
  for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );  
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticlesPostProcessing_Paraview_MPIIO_binary/types" ); 


  // Write field values to the binary buffer
  normU_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      normU = Norm( *(*particle)->getTranslationalVelocity() );
      nc = (*particle)->numberOfCells_PARAVIEW();
      for (i=0;i<nc;++i) write_double_binary( normU );
    }
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesPostProcessing_Paraview/NormU" ); 
  
  normOm_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      normOm = Norm( *(*particle)->getAngularVelocity() );
      nc = (*particle)->numberOfCells_PARAVIEW();
      for (i=0;i<nc;++i) write_double_binary( normOm );
    }
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesPostProcessing_Paraview/NormOm" ); 

  coord_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      coordNum = double((*particle)->getCoordinationNumber());
      nc = (*particle)->numberOfCells_PARAVIEW();
      for (i=0;i<nc;++i) write_double_binary( coordNum );
    }
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesPostProcessing_Paraview/CoordNumb" ); 
  total_offset = OFFSET;  
  
  int* total_binary_offset_per_proc = wrapper->AllGather_INT( total_offset );
  int* cumul_binary_offset_per_proc = new int[m_nprocs];
  cumul_binary_offset_per_proc[0] = 0;
  for (rank=1;rank<m_nprocs;rank++) 
    cumul_binary_offset_per_proc[rank] = cumul_binary_offset_per_proc[rank-1]
    	+ total_binary_offset_per_proc[rank-1];


  // Header per piece + general for 1st and last process
  ostringstream oss;
  if ( m_rank == 0 )
  {   
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
    oss << "compressor=\"vtkZLibDataCompressor\">\n";
    oss << "<UnstructuredGrid>\n";  
  }
  oss << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "offset=\"" << cumul_binary_offset_per_proc[m_rank] 
  	+ point_binary_offset << "\" format=\"appended\"></DataArray>\n"; 
  oss << "</Points>\n";
  oss << "<Cells>\n";
  oss << "<DataArray type=\"Int32\" Name=\"connectivity\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + connectivity_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"offsets\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + offsets_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"types\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + cellstype_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n";   
  oss << "</Cells>\n";
  oss << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">\n";
  oss << "<DataArray type=\"Float32\" Name=\"NormU\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + normU_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n";        
  oss << "<DataArray type=\"Float32\" Name=\"NormOm\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + normOm_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n"; 
  oss << "<DataArray type=\"Float32\" Name=\"CoordNumb\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + coord_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n"; 
  oss << "</CellData>\n";
  oss << "</Piece>\n";
  if ( m_rank == m_nprocs - 1 )
  { 
    oss << "</UnstructuredGrid>\n"; 
    oss << "<AppendedData encoding=\"raw\">\n" << "    _" ;  
  }            
  int header = int(oss.str().size());
  int* header_per_proc = wrapper->AllGather_INT( header );
  vector<int> mpifile_offsets( m_nprocs, 0 );
  mpifile_offsets[0] = 0;
  for (rank=1;rank<m_nprocs;rank++) 
    mpifile_offsets[rank] = mpifile_offsets[rank-1]
    	+ header_per_proc[rank-1] * sizeof(char);  
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss.str().c_str(), 
  	header, MPI_CHAR, &status ); 


  // Write the binary buffers
  int starting_point = mpifile_offsets[m_nprocs-1] 
  	+ header_per_proc[m_nprocs-1] * sizeof(char); 
  for (rank=0;rank<m_nprocs;rank++) 
    mpifile_offsets[rank] = starting_point
    	+ cumul_binary_offset_per_proc[rank] * sizeof(char);	
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], BUFFER, OFFSET, 
	MPI_CHAR, &status );	 	
  delete [] BUFFER ; BUFFER = 0 ;
  ALLOCATED = 0 ;
  OFFSET = 0 ;


  // Closing text
  ostringstream oss2;
  oss2 << "\n</AppendedData>\n";
  oss2 << "</VTKFile>\n";
  header = int(oss2.str().size());
  starting_point = mpifile_offsets[m_nprocs-1] 
  	+ total_binary_offset_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, starting_point, oss2.str().c_str(), header, 
    	MPI_CHAR, &status );

  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  delete [] total_binary_offset_per_proc;
  delete [] cumul_binary_offset_per_proc; 
  delete [] header_per_proc;    
}




// ----------------------------------------------------------------------------
// Writes spherical particles data in a vector form containing the
// center of mass coordinates of each particle in a single MPI file in text 
// mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeSpheresPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{

}




// ----------------------------------------------------------------------------
// Writes spherical particles data in a vector form containing the
// center of mass coordinates of each particle in a single MPI file in binary 
// mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeSpheresPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{

}




// ----------------------------------------------------------------------------
// Writes particle translational and angular velocity vectors in a single MPI 
// file in text mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticleVelocityVectorsPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{

}




// ----------------------------------------------------------------------------
// Writes particle translational and angular velocity vectors in a single MPI 
// file in binary mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticleVelocityVectorsPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{

}




// ----------------------------------------------------------------------------
// Writes contact force vectors in a single MPI file in text mode with MPI I/O 
// routines 
void ParaviewPostProcessingWriter::
	writeContactForceVectorsPostProcessing_Paraview_MPIIO_text(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites )
{

}




// ----------------------------------------------------------------------------
// Writes contact force vectors in a single MPI file in binary mode with MPI 
// I/O routines 
void ParaviewPostProcessingWriter::
	writeContactForceVectorsPostProcessing_Paraview_MPIIO_binary(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites )
{

}
