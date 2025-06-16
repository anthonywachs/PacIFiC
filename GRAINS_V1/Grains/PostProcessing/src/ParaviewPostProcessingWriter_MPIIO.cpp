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
	writeParticles_Paraview_MPIIO_text(
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
  char fmt[8] = "%12.5e ";
  char endfmt[8] = "%12.5e\n";
  const int charspernum = 13;
  size_t counter = 0;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  double nu, nom;     
  
  // Create a MPI datatype to write doublea as strings with a given format
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
  delete [] out_length_per_proc;  
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  delete [] max_offset_per_proc;


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
  delete [] out_length_per_proc;    
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;


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
  delete [] out_length_per_proc;
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
   
   
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
  delete [] out_length_per_proc;      
}




// ----------------------------------------------------------------------------
// Writes particles data in a single MPI file in binary mode with 
// MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticles_Paraview_MPIIO_binary(
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
  	"writeParticles_Paraview_MPIIO_binary/Points" );


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
  	"writeParticles_Paraview_MPIIO_binary/connectivity" );

  // Write offsets to the binary buffer  
  offsets_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, int(offsets.size()) ) ;  
  for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );  
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticles_Paraview_MPIIO_binary/offsets" ); 


  // Write cell types to the binary buffer  
  cellstype_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;  
  for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );  
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticles_Paraview_MPIIO_binary/types" ); 


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
  	"writeParticles_Paraview/NormU" ); 
  
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
  	"writeParticles_Paraview/NormOm" ); 

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
  	"writeParticles_Paraview/CoordNumb" ); 
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
	writeParticlesAsGlyph_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  MPI_Datatype num_as_string;
  int nbpts = 0, i, coordNum;
  list<Particle*>::const_iterator particle;
  char fmt[8] = "%12.5e ";
  char endfmt[8] = "%12.5e\n";
  const int charspernum = 13;
  size_t counter = 0;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  Point3 gc; 
  Vector3 vec;
  Quaternion qrot;   
  double nu, nom;       
  
  // Create a MPI datatype to write doublea as strings with a given format
  MPI_Type_contiguous( charspernum, MPI_CHAR, &num_as_string ); 
  MPI_Type_commit( &num_as_string );

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
      nbpts++;
  
  int total_nbpts = wrapper->sum_INT( nbpts ); 
  int* nbpts_per_proc = wrapper->AllGather_INT( nbpts );       


  // File header
  ostringstream oss;
  oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">\n";
  oss << "<UnstructuredGrid>\n";
  oss << "<Piece NumberOfPoints=\"" << total_nbpts << "\""
    	<< " NumberOfCells=\"" << ( total_nbpts ? "1" : "0" ) << "\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "format=\"ascii\">";
  oss << "\n";        
  int header = int(oss.str().size());
  if ( m_rank == 0 )
    MPI_File_write( file, oss.str().c_str(), header, MPI_CHAR, &status );

  // Write point coordinates to the MPI file
  char* coord = new char[ 3 * nbpts * charspernum + 1 ];
  list<Point3> ppp;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      gc = *(*particle)->getPosition();
      if ( PPTranslation ) gc += *PPTranslation ;
      sprintf( &coord[3*counter*charspernum], fmt, gc[X] );
      sprintf( &coord[(3*counter+1)*charspernum], fmt, gc[Y] );
      sprintf( &coord[(3*counter+2)*charspernum], endfmt, gc[Z] );
      ++counter;
    }    

  vector<int> mpifile_offsets( m_nprocs, 0 );
  mpifile_offsets[0] = header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 3 * nbpts, 
    	num_as_string, &status ); 
 

  // Connectivity, offsets and types and header for orientation
  ostringstream oss2;
  oss2 << "</DataArray>\n";
  oss2 << "</Points>\n";
  oss2 << "<Cells>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  oss2 << "0</DataArray>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
  oss2 << "1</DataArray>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">";
  oss2 << "1</DataArray>\n";
  oss2 << "</Cells>\n"; 
  oss2 << "<PointData ";
  oss2 << "Scalars=\"NormU,NormOm,CoordNumb\">\n";
  oss2 << "<DataArray Name=\"Quaternion\" "
    << "NumberOfComponents=\"4\" type=\"Float32\" format=\"ascii\">\n";         
  header = int(oss2.str().size());
  int mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss2.str().c_str(), header, 
    	MPI_CHAR, &status );
    
 
  // Quaternion
  delete [] coord;
  coord = new char[ 4 * nbpts * charspernum + 1 ];  
  counter = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      qrot = *(*particle)->getQuaternionRotation();
      sprintf( &coord[4*counter*charspernum], fmt, qrot[W] );      
      sprintf( &coord[(4*counter+1)*charspernum], fmt, qrot[X] );
      sprintf( &coord[(4*counter+2)*charspernum], fmt, qrot[Y] );
      sprintf( &coord[(4*counter+3)*charspernum], endfmt, qrot[Z] );
      ++counter;
    }    

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 4 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 4 * nbpts, 
    	num_as_string, &status );
	
  delete [] coord;


  // Norm of translational velocity
  ostringstream oss5;
  oss5 << "</DataArray>\n";
  oss5 << "<DataArray type=\"Float32\" Name=\"NormU\" format=\"ascii\">\n";    
  header = int(oss5.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 4 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss5.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  char* scalar = new char[ nbpts * charspernum + 1 ];
  counter = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nu = Norm( *(*particle)->getTranslationalVelocity() );
      sprintf( &scalar[counter*charspernum], fmt, nu );
      ++counter;	
    }

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * charspernum * sizeof(char) ;
 
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], scalar, nbpts, 
    	num_as_string, &status ); 
  
  
  // Norm of angular velocity
  ostringstream oss6;
  oss6 << "\n</DataArray>\n";  
  oss6 << "<DataArray type=\"Float32\" Name=\"NormOm\" format=\"ascii\">\n";    
  header = int(oss6.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss6.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  counter = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      nom = Norm( *(*particle)->getAngularVelocity() );
      sprintf( &scalar[counter*charspernum], fmt, nom );
      ++counter;	
    }

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * charspernum * sizeof(char) ;
 
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], scalar, nbpts, 
    	num_as_string, &status );
  
  
  // Coordination number
  ostringstream oss7;
  oss7 << "\n</DataArray>\n";  
  oss7 << "<DataArray type=\"Float32\" Name=\"CoordNumb\" ";        
  oss7 << "format=\"ascii\">\n";    
  header = int(oss7.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss7.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  ostringstream* oss_out = new ostringstream;
  for (particle=particles->begin();
    	particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      coordNum = (*particle)->getCoordinationNumber();
      *oss_out << coordNum << " ";
    }
  int out_length = int(oss_out->str().size());
  int* out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  delete [] scalar;
   
   
  // Closing text
  ostringstream oss8;
  oss8 << "\n</DataArray>\n";
  oss8 << "</PointData>\n";
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
  delete [] out_length_per_proc;     
}




// ----------------------------------------------------------------------------
// Writes spherical particles data in a vector form containing the
// center of mass coordinates of each particle in a single MPI file in binary 
// mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticlesAsGlyph_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  int nbpts = 0, rank;
  list<Particle*>::const_iterator particle;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  double normU, normOm, coordNum;
  int point_binary_offset = 0, orientation_binary_offset, 
  	normU_binary_offset, normOm_binary_offset, coord_binary_offset, 
	total_offset; 
  Point3 gc;
  Vector3 vec;
  Quaternion qrot;      	    

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
      nbpts++;


  // Write point coordinates to the binary buffer
  start_output_binary( sizeof_Float32, 3 * nbpts ) ;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      gc = *(*particle)->getPosition();
      if ( PPTranslation ) gc += *PPTranslation ;
      for (int comp=0;comp<3;++comp) write_double_binary( gc[comp] ) ;
    }   
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeParticlesAsGlyph_Paraview_MPIIO_binary/Points" );  


  // Write quaternion to the the binary buffer
  orientation_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, 4 * nbpts );  
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      qrot = *(*particle)->getQuaternionRotation();
      write_double_binary( qrot[W] ) ;
      for (int comp=0;comp<3;++comp) write_double_binary( qrot[comp] ) ;
    }    
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesAsGlyph_Paraview_MPIIO_binary/Orientation" ); 


  // Write field values to the binary buffer
  normU_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, nbpts );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      normU = Norm( *(*particle)->getTranslationalVelocity() );
      write_double_binary( normU );
    }
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesAsGlyph_Paraview_MPIIO_binary/NormU" ); 
  
  normOm_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, nbpts );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      normOm = Norm( *(*particle)->getAngularVelocity() );
      write_double_binary( normOm );
    }
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesAsGlyph_Paraview_MPIIO_binary/NormOm" ); 

  coord_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, nbpts );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      coordNum = double((*particle)->getCoordinationNumber());
      write_double_binary( coordNum );
    }
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeParticlesAsGlyph_Paraview_MPIIO_binary/CoordNumb" ); 
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
    	<< " NumberOfCells=\"" << ( nbpts ? "1" : "0" ) << "\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "offset=\"" << cumul_binary_offset_per_proc[m_rank] 
  	+ point_binary_offset << "\" format=\"appended\"></DataArray>\n"; 
  oss << "</Points>\n";
  oss << "<Cells>\n";
  oss << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  oss << "0</DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
  oss << "1</DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">";
  oss << "1</DataArray>\n";
  oss << "</Cells>\n"; 
  oss << "<PointData ";
  oss << "Scalars=\"NormU,NormOm,CoordNumb\">\n";
  oss << "<DataArray Name=\"Quaternion\" NumberOfComponents=\"4\" " <<
  	"type=\"Float32\" offset=\"" << cumul_binary_offset_per_proc[m_rank] 
	+ orientation_binary_offset << "\" format=\"appended\"></DataArray>\n";
  oss << "<DataArray type=\"Float32\" Name=\"NormU\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + normU_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n";        
  oss << "<DataArray type=\"Float32\" Name=\"NormOm\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + normOm_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n"; 
  oss << "<DataArray type=\"Float32\" Name=\"CoordNumb\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] + coord_binary_offset 
	<< "\" format=\"appended\"></DataArray>\n"; 
  oss << "</PointData>\n";
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
// Writes particle translational and angular velocity vectors in a single MPI 
// file in text mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticleVelocityVectors_Paraview_MPIIO_text(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  MPI_Datatype num_as_string;
  int nbpts = 0, i;
  list<Particle*>::const_iterator particle;
  char fmt[8] = "%12.5e ";
  char endfmt[8] = "%12.5e\n";
  const int charspernum = 13;
  size_t counter = 0;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  Point3 gc; 
  Vector3 const* vec;      
  
  // Create a MPI datatype to write doublea as strings with a given format
  MPI_Type_contiguous( charspernum, MPI_CHAR, &num_as_string ); 
  MPI_Type_commit( &num_as_string );

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
      nbpts++;
  
  int total_nbpts = wrapper->sum_INT( nbpts ); 
  int* nbpts_per_proc = wrapper->AllGather_INT( nbpts );       


  // File header
  ostringstream oss;
  oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">\n";
  oss << "<UnstructuredGrid>\n";
  oss << "<Piece NumberOfPoints=\"" << total_nbpts << "\""
    	<< " NumberOfCells=\"0\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "format=\"ascii\">";
  oss << "\n";        
  int header = int(oss.str().size());
  if ( m_rank == 0 )
    MPI_File_write( file, oss.str().c_str(), header, MPI_CHAR, &status );

  // Write point coordinates to the MPI file
  char* coord = new char[ 3 * nbpts * charspernum + 1 ];
  list<Point3> ppp;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      gc = *(*particle)->getPosition();
      if ( PPTranslation ) gc += *PPTranslation ;
      sprintf( &coord[3*counter*charspernum], fmt, gc[X] );
      sprintf( &coord[(3*counter+1)*charspernum], fmt, gc[Y] );
      sprintf( &coord[(3*counter+2)*charspernum], endfmt, gc[Z] );
      ++counter;
    }    

  vector<int> mpifile_offsets( m_nprocs, 0 );
  mpifile_offsets[0] = header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 3 * nbpts, 
    	num_as_string, &status ); 
  

  // Connectivity, offsets and types and header for orientation
  ostringstream oss2;
  oss2 << "</DataArray>\n";
  oss2 << "</Points>\n";
  oss2 << "<Cells>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  oss2 << "0 0 0</DataArray>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
  oss2 << "3</DataArray>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">";
  oss2 << "5</DataArray>\n";
  oss2 << "</Cells>\n"; 
  oss2 << "<PointData Vectors=\"U,Omega\">\n";
  oss2 << "<DataArray Name=\"U\" "
    << "NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";         
  header = int(oss2.str().size());
  int mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss2.str().c_str(), header, 
    	MPI_CHAR, &status );
    
  
  // Translational velocity vector
  counter = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      vec = (*particle)->getTranslationalVelocity();
      sprintf( &coord[3*counter*charspernum], fmt, (*vec)[X] );
      sprintf( &coord[(3*counter+1)*charspernum], fmt, (*vec)[Y] );
      sprintf( &coord[(3*counter+2)*charspernum], endfmt, (*vec)[Z] );
      ++counter;
    }    

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 3 * nbpts, 
    	num_as_string, &status );	
	

  // Angular velocity
  ostringstream oss5;
  oss5 << "</DataArray>\n";
  oss5 << "<DataArray Name=\"Omega\" "
    << "NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";    
  header = int(oss5.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss5.str().c_str(), header, 
    	MPI_CHAR, &status ); 

  counter = 0;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      vec = (*particle)->getAngularVelocity();
      sprintf( &coord[3*counter*charspernum], fmt, (*vec)[X] );
      sprintf( &coord[(3*counter+1)*charspernum], fmt, (*vec)[Y] );
      sprintf( &coord[(3*counter+2)*charspernum], endfmt, (*vec)[Z] );
      ++counter;
    }    

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 3 * nbpts, 
    	num_as_string, &status );
	
  delete [] coord;
   
   
  // Closing text
  ostringstream oss8;
  oss8 << "</DataArray>\n";
  oss8 << "</PointData>\n";
  oss8 << "</Piece>\n";
  oss8 << "</UnstructuredGrid>\n";
  oss8 << "</VTKFile>\n";
  header = int(oss8.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss8.str().c_str(), header, 
    	MPI_CHAR, &status );   


  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  MPI_Type_free( &num_as_string );
  delete [] nbpts_per_proc;
}




// ----------------------------------------------------------------------------
// Writes particle translational and angular velocity vectors in a single MPI 
// file in binary mode with MPI I/O routines
void ParaviewPostProcessingWriter::
	writeParticleVelocityVectors_Paraview_MPIIO_binary(
  	list<Particle*> const* particles, string const& partFilename,
	bool const& forceForAllTag,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  int nbpts = 0, rank;
  list<Particle*>::const_iterator particle;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  int point_binary_offset = 0, translational_binary_offset, 
  	angular_binary_offset, total_offset; 
  Point3 gc;
  Vector3 const* vec;     	    

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
       ( (*particle)->getTag() != 2 || forceForAllTag ) )
      nbpts++;


  // Write point coordinates to the binary buffer
  start_output_binary( sizeof_Float32, 3 * nbpts ) ;
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      gc = *(*particle)->getPosition();
      if ( PPTranslation ) gc += *PPTranslation ;
      for (int comp=0;comp<3;++comp) write_double_binary( gc[comp] ) ;
    }   
  compress_segment_binary( CURRENT_LENGTH,	
	"writePVelocityVectorsPostProcessing_Paraview_MPIIO_binary/Points" );  


  // Translational velocity vector
  translational_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, 3 * nbpts );  
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      vec = (*particle)->getTranslationalVelocity();
      for (int comp=0;comp<3;++comp) write_double_binary( (*vec)[comp] ) ;
    }    
  compress_segment_binary( CURRENT_LENGTH, 
  	"writePVelocityVectorsPostProcessing_Paraview_MPIIO_binary/U" ); 


  // Write field values to the binary buffer
  angular_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, 3 * nbpts );
  for (particle=particles->begin();particle!=particles->end();particle++)
    if ( (*particle)->getActivity() == COMPUTE && 
	( (*particle)->getTag() != 2 || forceForAllTag ) )
    {
      vec = (*particle)->getAngularVelocity();
      for (int comp=0;comp<3;++comp) write_double_binary( (*vec)[comp] ) ;
    }  
  compress_segment_binary( CURRENT_LENGTH, 
  	"writePVelocityVectorsPostProcessing_Paraview_MPIIO_binary/Omega" ); 
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
    	<< " NumberOfCells=\"0\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "offset=\"" << cumul_binary_offset_per_proc[m_rank] 
  	+ point_binary_offset << "\" format=\"appended\"></DataArray>\n"; 
  oss << "</Points>\n";
  oss << "<Cells>\n";
  oss << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  oss << "0 0 0</DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
  oss << "3</DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">";
  oss << "5</DataArray>\n";
  oss << "</Cells>\n"; 
  oss << "<PointData Vectors=\"U,Omega\">\n";
  oss << "<DataArray Name=\"U\" NumberOfComponents=\"3\" " <<
  	"type=\"Float32\" offset=\"" << cumul_binary_offset_per_proc[m_rank] 
	+ translational_binary_offset << "\" format=\"appended\">"
	<< "</DataArray>\n";
  oss << "<DataArray Name=\"Omega\" NumberOfComponents=\"3\" " <<
  	"type=\"Float32\" offset=\"" << cumul_binary_offset_per_proc[m_rank] 
	+ angular_binary_offset << "\" format=\"appended\">"
	<< "</DataArray>\n";
  oss << "</PointData>\n";
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
// Writes contact force vectors in a single MPI file in text mode with MPI I/O 
// routines 
void ParaviewPostProcessingWriter::
	writeContactForceVectors_Paraview_MPIIO_text(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  MPI_Datatype num_as_string;
  int nbpts = 0, i;
  char fmt[8] = "%12.5e ";
  char endfmt[8] = "%12.5e\n";
  const int charspernum = 13;
  size_t counter = 0;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  vector<struct PointForcePostProcessing> const* pallContacts = 
  	LC->getPPForces();
  Point3 pt;
        
  
  // Create a MPI datatype to write doublea as strings with a given format
  MPI_Type_contiguous( charspernum, MPI_CHAR, &num_as_string ); 
  MPI_Type_commit( &num_as_string );

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  nbpts = int(LC->getNbPPForces());
  
  int total_nbpts = wrapper->sum_INT( nbpts ); 
  int* nbpts_per_proc = wrapper->AllGather_INT( nbpts );       


  // File header
  ostringstream oss;
  oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">\n";
  oss << "<UnstructuredGrid>\n";
  oss << "<Piece NumberOfPoints=\"" << total_nbpts << "\""
    	<< " NumberOfCells=\"0\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "format=\"ascii\">\n";        
  int header = int(oss.str().size());
  if ( m_rank == 0 )
    MPI_File_write( file, oss.str().c_str(), header, MPI_CHAR, &status );


  // Write point coordinates to the MPI file
  char* coord = new char[ 3 * nbpts * charspernum + 1 ];
  for (i=0;i<nbpts;++i)
  {
    pt = (*pallContacts)[i].geometricPointOfContact;
    if ( PPTranslation ) pt += *PPTranslation ;
    sprintf( &coord[3*counter*charspernum], fmt, pt[X] );
    sprintf( &coord[(3*counter+1)*charspernum], fmt, pt[Y] );
    sprintf( &coord[(3*counter+2)*charspernum], endfmt, pt[Z] );
    ++counter;
  }         

  vector<int> mpifile_offsets( m_nprocs, 0 );
  mpifile_offsets[0] = header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 3 * nbpts, 
    	num_as_string, &status ); 
  

  // Connectivity, offsets and types and header for orientation
  ostringstream oss2;
  oss2 << "</DataArray>\n";
  oss2 << "</Points>\n";
  oss2 << "<Cells>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  oss2 << "0 0 0</DataArray>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
  oss2 << "3</DataArray>\n";
  oss2 << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">";
  oss2 << "5</DataArray>\n";
  oss2 << "</Cells>\n"; 
  oss2 << "<PointData Vectors=\"Force\">\n";
  oss2 << "<DataArray Name=\"Force\" "
    << "NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";         
  header = int(oss2.str().size());
  int mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss2.str().c_str(), header, 
    	MPI_CHAR, &status );
    
  
  // Force vector
  counter = 0;
  for (i=0;i<nbpts;++i)
  {
    sprintf( &coord[3*counter*charspernum], fmt, 
    	(*pallContacts)[i].contactForceComp0[X] );
    sprintf( &coord[(3*counter+1)*charspernum], fmt, 
    	(*pallContacts)[i].contactForceComp0[Y] );
    sprintf( &coord[(3*counter+2)*charspernum], endfmt, 
    	(*pallContacts)[i].contactForceComp0[Z] );
    ++counter;
  }  

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], coord, 3 * nbpts, 
    	num_as_string, &status );	
	
  delete [] coord;
   
   
  // Closing text
  ostringstream oss8;
  oss8 << "</DataArray>\n";
  oss8 << "</PointData>\n";
  oss8 << "</Piece>\n";
  oss8 << "</UnstructuredGrid>\n";
  oss8 << "</VTKFile>\n";
  header = int(oss8.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss8.str().c_str(), header, 
    	MPI_CHAR, &status );   


  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  MPI_Type_free( &num_as_string );
  delete [] nbpts_per_proc;
}




// ----------------------------------------------------------------------------
// Writes contact force vectors in a single MPI file in binary mode with MPI 
// I/O routines 
void ParaviewPostProcessingWriter::
	writeContactForceVectors_Paraview_MPIIO_binary(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& partFilename, double const& time,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  int nbpts = 0, i, rank, point_binary_offset = 0, force_binary_offset, 
  	total_offset; 
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  vector<struct PointForcePostProcessing> const* pallContacts = 
  	LC->getPPForces();
  Point3 pt;    	    

  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + partFilename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  nbpts = int(LC->getNbPPForces());


  // Write point coordinates to the binary buffer
  start_output_binary( sizeof_Float32, 3 * nbpts );    
  for (i=0;i<nbpts;++i)
  {
    pt = (*pallContacts)[i].geometricPointOfContact;
    if ( PPTranslation ) pt += *PPTranslation ;
    for (int comp=0;comp<3;++comp) write_double_binary( pt[comp] ) ;
  }
  compress_segment_binary( CURRENT_LENGTH,	
	"writeContactForceVectors_Paraview_MPIIO_binary/Points" );


  // Contact force vector
  force_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, 3 * nbpts );  
  for (i=0;i<nbpts;++i)
     for (int comp=0;comp<3;++comp) write_double_binary( 
     	(*pallContacts)[i].contactForceComp0[comp] ) ;    
  compress_segment_binary( CURRENT_LENGTH, 
  	"writeContactForceVectors_Paraview_MPIIO_binary/Force" ); 
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
    	<< " NumberOfCells=\"0\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "offset=\"" << cumul_binary_offset_per_proc[m_rank] 
  	+ point_binary_offset << "\" format=\"appended\"></DataArray>\n"; 
  oss << "</Points>\n";
  oss << "<Cells>\n";
  oss << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  oss << "0 0 0</DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
  oss << "3</DataArray>\n";
  oss << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">";
  oss << "5</DataArray>\n";
  oss << "</Cells>\n"; 
  oss << "<PointData Vectors=\"Force\">\n";
  oss << "<DataArray Name=\"Force\" NumberOfComponents=\"3\" " <<
  	"type=\"Float32\" offset=\"" << cumul_binary_offset_per_proc[m_rank] 
	+ force_binary_offset << "\" format=\"appended\">"
	<< "</DataArray>\n";
  oss << "</PointData>\n";
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
// Writes contact force chain network in a single MPI file in text mode with 
// MPI I/O routines 
void ParaviewPostProcessingWriter::
	writeContactForceChains_Paraview_MPIIO_text(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& filename, double const& time,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  MPI_Datatype num_as_string;
  char fmt[8] = "%12.5e ";
  char endfmt[8] = "%12.5e\n";
  const int charspernum = 13;
  size_t counter = 0;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  vector<struct PointForcePostProcessing> const* pallContacts = 
  	LC->getPPForces();
  int i, j, nPPF = int(LC->getNbPPForces()), nbpts, mpifile_offset, rank;
  Point3 pt0, pt1;
  double norm = 0.;
  
  
  // Create a MPI datatype to write doublea as strings with a given format
  MPI_Type_contiguous( charspernum, MPI_CHAR, &num_as_string ); 
  MPI_Type_commit( &num_as_string );
  
  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + filename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file ); 
	

  // Numbers of Paraview points
  nbpts = 2 * nPPF;
  
  int total_nbpts = wrapper->sum_INT( nbpts ); 
  int* nbpts_per_proc = wrapper->AllGather_INT( nbpts ); 


  // File header
  ostringstream oss;
  oss << "<VTKFile type=\"PolyData\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\">\n";
  oss << "<PolyData>\n";
  oss << "<Piece NumberOfPoints=\"" << total_nbpts 
    << "\" NumberOfVerts=\"" << 0
    << "\" NumberOfLines=\"" << total_nbpts / 2
    << "\" NumberOfStrips=\"" << 0
    << "\" NumberOfPolys=\"" << 0 << "\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "format=\"ascii\">\n";        
  int header = int(oss.str().size());
  if ( m_rank == 0 )
    MPI_File_write( file, oss.str().c_str(), header, MPI_CHAR, &status );
    
  // Write point coordinates to the MPI file
  char* pts_coord = new char[ 3 * nbpts * charspernum + 1 ];
  for (i=0;i<nPPF;++i,counter++)
  {
    // Point 0
    pt0 = (*pallContacts)[i].PPptComp0;
    if ( PPTranslation ) pt0 += *PPTranslation ;
    sprintf( &pts_coord[6*counter*charspernum], fmt, pt0[X] );
    sprintf( &pts_coord[(6*counter+1)*charspernum], fmt, pt0[Y] );
    sprintf( &pts_coord[(6*counter+2)*charspernum], endfmt, pt0[Z] );     

    // Point 1
    pt1 = (*pallContacts)[i].PPptComp1;
    if ( PPTranslation ) pt1 += *PPTranslation ; 
    sprintf( &pts_coord[(6*counter+3)*charspernum], fmt, pt1[X] );
    sprintf( &pts_coord[(6*counter+4)*charspernum], fmt, pt1[Y] );
    sprintf( &pts_coord[(6*counter+5)*charspernum], endfmt, pt1[Z] );
  }

  vector<int> mpifile_offsets( m_nprocs, 0 );
  mpifile_offsets[0] = header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * 3 * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], pts_coord, 3 * nbpts, 
    	num_as_string, &status ); 

  delete [] pts_coord; 

  
  // Force magnitude
  ostringstream oss2;
  oss2 << "</DataArray>\n";
  oss2 << "</Points>\n";
  oss2 << "<PointData Scalars=\"ForceMag\">\n";
  oss2 << "<DataArray type=\"Float32\" Name=\"ForceMag\" ";        
  oss2 << "format=\"ascii\">\n";    
  header = int(oss2.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * 3 * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss2.str().c_str(), header, 
    	MPI_CHAR, &status ); 
	
  char* normF = new char[ nbpts * charspernum + 1 ];
  counter = 0;
  for (i=0;i<nPPF;++i)
  {
    norm = Norm( (*pallContacts)[i].contactForceComp0 );
    sprintf( &normF[counter*charspernum], fmt, norm );
    sprintf( &normF[(counter+1)*charspernum], fmt, norm );    
    counter += 2;    
  }

  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * charspernum * sizeof(char) ;
 
  MPI_File_write_at_all( file, mpifile_offsets[m_rank], normF, nbpts, 
    	num_as_string, &status ); 

  delete [] normF; 


  // Line connectivity
  ostringstream oss3;
  oss3 << "\n</DataArray>\n";
  oss3 << "</PointData>\n";
  oss3 << "<Lines>\n";
  oss3 << "<DataArray type=\"Int32\" Name=\"connectivity\" ";        
  oss3 << "format=\"ascii\">\n";    
  header = int(oss3.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ nbpts_per_proc[m_nprocs-1] * charspernum * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss3.str().c_str(), header, 
    	MPI_CHAR, &status );
	
  // Compute the point number shifts 
  vector<int> shift( m_nprocs, 0 );
  for (rank=1;rank<m_nprocs;rank++)
    for (j=0;j<rank;j++)
      shift[rank] += nbpts_per_proc[j];
      
  // Write connectivity to the MPI file with the point number shifts
  ostringstream* oss_out = new ostringstream;
  for (j=0;j<nbpts;j++)
    *oss_out << j + shift[m_rank] << " ";
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


  // Offsets
  ostringstream oss4;
  oss4 << "</DataArray>\n";
  oss4 << "<DataArray type=\"Int32\" Name=\"offsets\" ";        
  oss4 << "format=\"ascii\">\n";    
  header = int(oss4.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ out_length_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss4.str().c_str(), header, 
    	MPI_CHAR, &status );
	
  // Write offsets to the MPI file with the offset shifts
  oss_out = new ostringstream;
  for (j=0;j<nPPF;j++)
    *oss_out << 2 * j + shift[m_rank] << " ";
  if ( m_rank == m_nprocs - 1 ) *oss_out << "\n";  
  out_length = int(oss_out->str().size());
  delete [] out_length_per_proc;
  out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = mpifile_offset + header * sizeof(char);
  for (i=1;i<m_nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[m_rank], oss_out->str().c_str(), 
  	out_length, MPI_CHAR, &status );   

  delete oss_out;
  
  
  // Closing text
  ostringstream oss5;
  oss5 << "</DataArray>\n";
  oss5 << "</Lines>\n";
  oss5 << "</Piece>\n";
  oss5 << "</PolyData>\n";
  oss5 << "</VTKFile>\n";
  header = int(oss5.str().size());
  mpifile_offset = mpifile_offsets[m_nprocs-1] 
  	+ out_length_per_proc[m_nprocs-1] * sizeof(char);
  if ( m_rank == 0 )
    MPI_File_write_at( file, mpifile_offset, oss5.str().c_str(), header, 
    	MPI_CHAR, &status );  
                
	
  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  MPI_Type_free( &num_as_string );
  delete [] nbpts_per_proc;
  delete [] out_length_per_proc;    	    
}




// ----------------------------------------------------------------------------
// Writes contact force chain network in a single MPI file in binary mode with 
// MPI I/O routines 
void ParaviewPostProcessingWriter::
	writeContactForceChains_Paraview_MPIIO_binary(
  	list<Particle*> const* particles,
  	LinkedCell const* LC, string const& filename, double const& time,
	bool const& processwrites )
{
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status;
  Vector3 const* PPTranslation = 
  	GrainsExec::m_translationParaviewPostProcessing ;
  vector<struct PointForcePostProcessing> const* pallContacts = 
  	LC->getPPForces();
  int i, j, rank, point_binary_offset = 0, force_binary_offset, 
  	connectivity_binary_offset, offsets_binary_offset, total_offset, 
	nPPF = int(LC->getNbPPForces()), nbpts, comp; 	
  Point3 pt0, pt1;
  double norm = 0.;
  
  // Open the file 
  MPI_File_open( MPI_COMM_activeProc, ( m_ParaviewFilename_dir 
  	+ "/" + filename ).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  // Numbers of Paraview points
  nbpts = 2 * nPPF;
  
  int* nbpts_per_proc = wrapper->AllGather_INT( nbpts );
  
  
  // Write point coordinates to the binary buffer
  start_output_binary( sizeof_Float32, 3 * nbpts );    
  for (i=0;i<nPPF;++i)
  {
    // Point 0
    pt0 = (*pallContacts)[i].PPptComp0;
    if ( PPTranslation ) pt0 += *PPTranslation ;
    for (comp=0;comp<3;++comp) write_double_binary( pt0[comp] ) ; 

    // Point 1
    pt1 = (*pallContacts)[i].PPptComp1;
    if ( PPTranslation ) pt1 += *PPTranslation ; 
    for (comp=0;comp<3;++comp) write_double_binary( pt1[comp] ) ; 
  }
  compress_segment_binary( CURRENT_LENGTH,	
	"writeContactForceChains_Paraview_MPIIO_binary/Points" ); 


  // Compute the point number shifts 
  vector<int> shift( m_nprocs, 0 );
  for (rank=1;rank<m_nprocs;rank++)
    for (j=0;j<rank;j++)
      shift[rank] += nbpts_per_proc[j];

	
  // Line connectivity
  connectivity_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, nbpts ) ;  
  for (j=0;j<nbpts;j++) write_int_binary( j );
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeContactForceChains_Paraview_MPIIO_binary/connectivity" );
	

  // Offsets
  offsets_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, nPPF ) ; 
  for (j=0;j<nPPF;j++) write_int_binary( 2 * j );
  compress_segment_binary( CURRENT_LENGTH,	
  	"writeContactForceChains_Paraview_MPIIO_binary/offsets" );


  // Force magnitude
  force_binary_offset = OFFSET;
  start_output_binary( sizeof_Float32, nbpts );    
  for (i=0;i<nPPF;++i)
  {
    norm = Norm( (*pallContacts)[i].contactForceComp0 );
    write_double_binary( norm ) ;
    write_double_binary( norm ) ;       
  }  
  compress_segment_binary( CURRENT_LENGTH,	
	"writeContactForceChains_Paraview_MPIIO_binary/Force" );  
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
    oss << "<VTKFile type=\"PolyData\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
    oss << "compressor=\"vtkZLibDataCompressor\">\n";
    oss << "<PolyData>\n";  
  }
  oss << "<Piece NumberOfPoints=\"" << nbpts 
    << "\" NumberOfVerts=\"" << 0
    << "\" NumberOfLines=\"" << nPPF
    << "\" NumberOfStrips=\"" << 0
    << "\" NumberOfPolys=\"" << 0 << "\">\n";
  oss << "<Points>\n";
  oss << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  oss << "offset=\"" << cumul_binary_offset_per_proc[m_rank] 
  	+ point_binary_offset << "\" format=\"appended\"></DataArray>\n"; 
  oss << "</Points>\n";
  oss << "<Lines>\n";
  oss << "<DataArray Name=\"connectivity\" type=\"Int32\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] 
	+ connectivity_binary_offset << "\" format=\"appended\">"
	<< "</DataArray>\n";
  oss << "<DataArray Name=\"offsets\" type=\"Int32\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] 
	+ offsets_binary_offset << "\" format=\"appended\">"
	<< "</DataArray>\n";
  oss << "</Lines>\n";  
  oss << "<PointData Scalars=\"ForceMag\">\n";
  oss << "<DataArray Name=\"ForceMag\" type=\"Float32\" offset=\"" << 
  	cumul_binary_offset_per_proc[m_rank] 
	+ force_binary_offset << "\" format=\"appended\">"
	<< "</DataArray>\n";
  oss << "</PointData>\n";		    
  oss << "</Piece>\n";
  if ( m_rank == m_nprocs - 1 )
  { 
    oss << "</PolyData>\n"; 
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
  delete [] nbpts_per_proc;
  delete [] total_binary_offset_per_proc;
  delete [] cumul_binary_offset_per_proc; 
  delete [] header_per_proc;	  	         
}
