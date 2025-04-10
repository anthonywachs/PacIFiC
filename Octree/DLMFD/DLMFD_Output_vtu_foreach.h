/** 
# Paraview functions 
*/
#include <zlib.h>


# if dimension == 2
#   define NVERTCELL 4
#   define CELLTYPE 9 
# endif
# if dimension == 3
#   define NVERTCELL 8
#   define CELLTYPE 12  
# endif

/**
# init_percelldir
This function returns an array of 4 or 8 vectors that enables to compute the 
coordinates of the vertices of a cell.
The vertex coordinates are : cell.center.x + percelldir[j].x * Delta 
*/
coord* init_percelldir()
{ 
  # if dimension == 2
    coord* percelldir = (coord*) malloc( 4 * sizeof(coord) );
    // Vertex [] = [0,0]
    percelldir[0].x = - 0.5;
    percelldir[0].y = - 0.5;
    percelldir[0].z = 0.;     
    // Vertex [1,0]
    percelldir[1].x = 0.5;
    percelldir[1].y = - 0.5;
    percelldir[1].z = 0.;         
    // Vertex [1,1]
    percelldir[2].x = 0.5;
    percelldir[2].y = 0.5; 
    percelldir[2].z = 0.;        
    // Vertex [0,1]
    percelldir[3].x = - 0.5;
    percelldir[3].y = 0.5;
    percelldir[3].z = 0.;                
# endif
# if dimension == 3
    coord* percelldir = (coord*) malloc( 8 * sizeof(coord) );
    // Vertex [] = [0,0,0]
    percelldir[0].x = - 0.5;
    percelldir[0].y = - 0.5;
    percelldir[0].z = - 0.5;     
    // Vertex [1,0,0]
    percelldir[1].x = 0.5;
    percelldir[1].y = - 0.5;
    percelldir[1].z = - 0.5;         
    // Vertex [1,1,0]
    percelldir[2].x = 0.5;
    percelldir[2].y = 0.5; 
    percelldir[2].z = - 0.5;        
    // Vertex [0,1,0]
    percelldir[3].x = - 0.5;
    percelldir[3].y = 0.5;
    percelldir[3].z = - 0.5;      
    // Vertex [] = [0,0,1]
    percelldir[4].x = - 0.5;
    percelldir[4].y = - 0.5;
    percelldir[4].z = 0.5;     
    // Vertex [1,0,1]
    percelldir[5].x = 0.5;
    percelldir[5].y = - 0.5;
    percelldir[5].z = 0.5;         
    // Vertex [1,1,1]
    percelldir[6].x = 0.5;
    percelldir[6].y = 0.5; 
    percelldir[6].z = 0.5;        
    // Vertex [0,1,1]
    percelldir[7].x = - 0.5;
    percelldir[7].y = 0.5;
    percelldir[7].z = 0.5;          
# endif

  return percelldir;
}


char* BUFFER ;
int ALLOCATED ;
int OFFSET ;
int CURRENT_LENGTH ;
int sizeof_Float32 = 4 ;
int sizeof_Int32 = 4 ;



/**
# output_pvtu_ascii
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_ascii_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_ascii( scalar* list, vector* vlist, FILE* fp, 
	char* subname )
{
  fputs( "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp );
  fputs( "<PUnstructuredGrid GhostLevel=\"0\">\n", fp );
  fputs( "<PCellData Scalars=\"scalars\">\n", fp );
  for (scalar s in list) 
  {
    fprintf( fp, "<PDataArray type=\"%s\" "
    	"Name=\"%s\" format=\"ascii\">\n", PARAVIEW_DATANAME, s.name );
    fputs( "</PDataArray>\n", fp );
  }
  for (vector v in vlist) 
  {
    fprintf( fp, "<PDataArray type=\"%s\" "
    	"NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", 
	PARAVIEW_DATANAME, v.x.name );
    fputs( "</PDataArray>\n", fp );
  }
  fputs( "</PCellData>\n", fp );
  fputs( "<PPoints>\n", fp );
  fprintf( fp, "<PDataArray type=\"%s\" "
  	"NumberOfComponents=\"3\" format=\"ascii\">\n", PARAVIEW_DATANAME );
  fputs( "</PDataArray>\n", fp );
  fputs( "</PPoints>\n", fp );

  for (unsigned int i = 0; i < npe(); i++)
    fprintf( fp, "<Piece Source=\"%s_%d.vtu\"/> \n", subname, i );

  fputs( "</PUnstructuredGrid>\n", fp );
  fputs( "</VTKFile>\n", fp );
}




/**
# output_vtu_asci_foreach
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one 
writes one *.vtu file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel.  
*/
void output_vtu_ascii_foreach( scalar* list, vector* vlist, 
	FILE* fp, bool linear )
{
# if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads( 1 );
# endif

  vertex scalar marker[];
  // We use uint32_t for the 3 following variables because the number
  // of cells/vertices on a single proc never exceeds the limit of 4,294,967,295
  uint32_t no_points = 0, no_cells = 0, vertexnum = 0 ;

  // In case of periodicity
  scalar per_mask[];
  foreach(serial, noauto) per_mask[] = 1.; 
  coord* percelldir = NULL;
  if ( Period.x || Period.y || Period.z )
  {  
    percelldir = init_percelldir();
    foreach(serial, noauto)
    {
      if ( Period.x )
        if ( x + Delta > X0 + L0 || x - Delta < X0 )
          per_mask[] = 0.;
#     if dimension > 1
        if ( Period.y )
          if ( y + Delta > Y0 + L0 ||  y - Delta < Y0 )
            per_mask[] = 0.;
#     endif
#     if dimension > 2
        if ( Period.z )
          if ( z + Delta > Z0 + L0 || z - Delta < Z0 )
            per_mask[] = 0.;
#     endif
    }
  }

  foreach_vertex(serial, noauto)
  {
    marker[] = vertexnum++;
    no_points += 1;
  }
  foreach(serial, noauto)
  {
    // Additional (duplicated) vertices per cell that have a face belonging to a
    // periodic face of the whole domain
    if ( per_mask[] < 0.5 )
      no_points += NVERTCELL;
    
    // Number of cells on this process
    no_cells += 1;
  }

  fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
  	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp );
  fputs( "<UnstructuredGrid>\n", fp );
  fprintf( fp, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n", 
  	no_points, no_cells );
  fputs( "<CellData Scalars=\"scalars\">\n", fp );
  for (scalar s in list) 
  {
    fprintf( fp, "<DataArray type=\"%s\" Name=\"%s\" "
    	"format=\"ascii\">\n", PARAVIEW_DATANAME, s.name );
    foreach(serial, noauto)
      fprintf( fp, "%g\n", val(s) );
    fputs( "</DataArray>\n", fp );
  }
  for (vector v in vlist) 
  {
    fprintf( fp, "<DataArray type=\"%s\" "
    	"NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", 
	PARAVIEW_DATANAME, v.x.name );
    foreach(serial, noauto)
    {
#     if dimension == 2
        fprintf( fp, "%g %g 0.\n", val(v.x), val(v.y) );
#     endif
#     if dimension == 3
        fprintf( fp, "%g %g %g\n", val(v.x), val(v.y), val(v.z) );
#     endif
    }
    fputs( "</DataArray>\n", fp );
  }
  fputs( "</CellData>\n", fp );
  fputs( "<Points>\n", fp );
  fprintf( fp, "<DataArray type=\"%s\" "
  	"NumberOfComponents=\"3\" format=\"ascii\">\n", PARAVIEW_DATANAME );
  foreach_vertex(serial, noauto)
  {
#   if dimension == 2
      fprintf( fp, "%g %g 0\n", x, y );
#   endif
#   if dimension == 3
      fprintf( fp, "%g %g %g\n", x, y, z );
#   endif
  }
  // Additional duplicated vertices in case of periodicity
  coord supvertex = {0., 0., 0.};
  if ( Period.x || Period.y || Period.z )
  {
    foreach(serial, noauto)
      if ( per_mask[] < 0.5 )
      {
#       if dimension == 2
          for (size_t j=0;j<4;j++)
	  {
	    supvertex.x = x + percelldir[j].x * Delta;
	    supvertex.y = y + percelldir[j].y * Delta;	    
            fprintf( fp, "%g %g 0.\n", supvertex.x, supvertex.y );	  
	  }		
#       endif
#       if dimension == 3
          for (size_t j=0;j<8;j++)
	  {
	    supvertex.x = x + percelldir[j].x * Delta;
	    supvertex.y = y + percelldir[j].y * Delta;	
	    supvertex.z = z + percelldir[j].z * Delta;	    	    
            fprintf( fp, "%g %g %g\n", supvertex.x, supvertex.y, supvertex.z );
	  }	
#       endif    
      }
    free( percelldir );       
  }  
  fputs( "</DataArray>\n", fp );
  fputs( "</Points>\n", fp );
  fputs( "<Cells>\n", fp );
  fputs( "<DataArray type=\"UInt32\" Name=\"connectivity\" "
  	"format=\"ascii\">\n", fp );
  foreach(serial, noauto)
    if ( per_mask[] )
    {
#     if dimension == 2
        uint32_t ape1 = (uint32_t)marker[];
        uint32_t ape2 = (uint32_t)marker[1,0];
        uint32_t ape3 = (uint32_t)marker[1,1];
        uint32_t ape4 = (uint32_t)marker[0,1];
        fprintf( fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4 );
#     endif
#     if dimension == 3
        uint32_t ape1 = (uint32_t)marker[];
        uint32_t ape2 = (uint32_t)marker[1,0,0];
        uint32_t ape3 = (uint32_t)marker[1,1,0];
        uint32_t ape4 = (uint32_t)marker[0,1,0];
        uint32_t ape5 = (uint32_t)marker[0,0,1];
        uint32_t ape6 = (uint32_t)marker[1,0,1];
        uint32_t ape7 = (uint32_t)marker[1,1,1];
        uint32_t ape8 = (uint32_t)marker[0,1,1];
        fprintf( fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, 
      		ape6, ape7, ape8 );
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        uint32_t ape1 = vertexnum++;
        uint32_t ape2 = vertexnum++;
        uint32_t ape3 = vertexnum++;
        uint32_t ape4 = vertexnum++;
        fprintf( fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4 );
#     endif
#     if dimension == 3
        uint32_t ape1 = vertexnum++;
        uint32_t ape2 = vertexnum++;
        uint32_t ape3 = vertexnum++;
        uint32_t ape4 = vertexnum++;
        uint32_t ape5 = vertexnum++;
        uint32_t ape6 = vertexnum++;
        uint32_t ape7 = vertexnum++;
        uint32_t ape8 = vertexnum++;
        fprintf( fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, 
      		ape6, ape7, ape8 );
#     endif
    }
  fputs( "</DataArray>\n", fp );
  fputs( "<DataArray type=\"UInt32\" Name=\"offsets\" "
  	"format=\"ascii\">\n", fp );
  for (uint32_t i = 1; i < no_cells+1; i++)
    fprintf( fp, "%u ", i * NVERTCELL );
  fprintf( fp, "\n" );    
  fputs( "</DataArray>\n", fp );
  fputs( "<DataArray type=\"UInt8\" Name=\"types\" "
  	"format=\"ascii\">\n", fp );
  foreach(serial, noauto)
    fprintf( fp, "%u ", CELLTYPE );
  fprintf( fp, "\n" );    
  fputs( "</DataArray>\n", fp );
  fputs( "</Cells>\n", fp );
  fputs( "</Piece>\n", fp );
  fputs( "</UnstructuredGrid>\n", fp );
  fputs( "</VTKFile>\n", fp );
  fflush( fp );
# if defined(_OPENMP)
    omp_set_num_threads( num_omp );
#  endif
}




/**
# output_pvtu_bin
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_bin_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_bin( scalar* list, vector* vlist,  FILE* fp, char* subname )
{
  fputs( "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp );
  fputs( " <PUnstructuredGrid GhostLevel=\"0\">\n", fp );
  fputs( "<PCellData Scalars=\"scalars\">\n", fp );
  for (scalar s in list) 
  {
    fprintf( fp,"<PDataArray type=\"%s\" Name=\"%s\" "
    	"format=\"appended\">\n", PARAVIEW_DATANAME, s.name );
    fputs( "</PDataArray>\n", fp );
  }
  for (vector v in vlist) 
  {
    fprintf( fp,"<PDataArray type=\"%s\" NumberOfComponents=\"3\""
    	" Name=\"%s\" format=\"appended\">\n", PARAVIEW_DATANAME, v.x.name );
    fputs( "</PDataArray>\n", fp );
  }
  fputs( "</PCellData>\n", fp );
  fputs( "<PPoints>\n", fp );
  fprintf( fp, "<PDataArray type=\"%s\" NumberOfComponents=\"3\" "
  	"format=\"ascii\">\n", PARAVIEW_DATANAME );
  fputs( " </PDataArray>\n", fp );
  fputs( "</PPoints>\n", fp );

  for (int i = 0; i < npe(); i++)
    fprintf( fp, "<Piece Source=\"%s_%d.vtu\"/> \n", subname, i );

  fputs( "</PUnstructuredGrid>\n", fp );
  fputs( "</VTKFile>\n", fp );
}




/**
# output_vtu_bin_foreach
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on binary format. If one 
writes one *.vtu file per PID process this function may be combined with
output_pvtu_bin() above to read in parallel. 
*/
void output_vtu_bin_foreach( scalar* list, vector* vlist, FILE* fp, 
	bool linear )
{
# if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads( 1 );
# endif

  vertex scalar marker[];
  // We use uint32_t for the 3 following variables because the number
  // of cells/vertices on a single proc never exceeds the limit of 4,294,967,295 
  uint32_t no_points = 0, no_cells = 0, vertexnum = 0 ;

  // In case of periodicity
  scalar per_mask[];
  foreach(serial, noauto) per_mask[] = 1.; 
  coord* percelldir = NULL;
  if ( Period.x || Period.y || Period.z )
  {  
    percelldir = init_percelldir();
    foreach(serial, noauto)
    {
      if ( Period.x )
        if ( x + Delta > X0 + L0 || x - Delta < X0 )
          per_mask[] = 0.;
#     if dimension > 1
        if ( Period.y )
          if ( y + Delta > Y0 + L0 ||  y - Delta < Y0 )
            per_mask[] = 0.;
#     endif
#     if dimension > 2
        if ( Period.z )
          if ( z + Delta > Z0 + L0 || z - Delta < Z0 )
            per_mask[] = 0.;
#     endif
    }
  }

  foreach_vertex(serial, noauto)
  {
    marker[] = vertexnum++;
    no_points += 1;
  }
  foreach(serial, noauto)
  {
    // Additional (duplicated) vertices per cell that have a face belonging to a
    // periodic face of the whole domain
    if ( per_mask[] < 0.5 ) 
      no_points += NVERTCELL;
    
    // Number of cells on this process
    no_cells += 1;
  }

  fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp );
  fputs( "<UnstructuredGrid>\n", fp );
  fprintf( fp,"<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n", 
  	no_points, no_cells );
  fputs( "<CellData Scalars=\"scalars\">\n", fp );

  uint64_t count = 0;
  for (scalar s in list) 
  {
    fprintf( fp,"<DataArray type=\"%s\" Name=\"%s\" "
    	"format=\"appended\" offset=\"%lu\">\n", PARAVIEW_DATANAME, 
	s.name, count );
    count += (uint64_t)no_cells * sizeof(PARAVIEW_DATATYPE) + sizeof(uint64_t);
    fputs( "</DataArray>\n", fp );
  }
  for (vector v in vlist) 
  {
    fprintf( fp,"<DataArray type=\"%s\" Name=\"%s\" "
    	"NumberOfComponents=\"3\" format=\"appended\" offset=\"%lu\">\n", 
	PARAVIEW_DATANAME, v.x.name, count );
    count += 3 * (uint64_t)no_cells * sizeof(PARAVIEW_DATATYPE) 
    	+ sizeof(uint64_t);
    fputs( "</DataArray>\n", fp );
  }
  fputs( "</CellData>\n", fp );
  fputs( "<Points>\n", fp );
  fprintf( fp, "<DataArray type=\"%s\" NumberOfComponents=\"3\" "
  	"format=\"appended\" offset=\"%lu\">\n", PARAVIEW_DATANAME, count );
  count += 3 * (uint64_t)no_points * sizeof(PARAVIEW_DATATYPE) 
  	+ sizeof(uint64_t);
  fputs( "</DataArray>\n", fp );
  fputs( "</Points>\n", fp );
  fputs( "<Cells>\n", fp );
  fprintf( fp,"<DataArray type=\"UInt32\" Name=\"connectivity\" "
  	"format=\"appended\" offset=\"%lu\"/>\n", count );
  count += (uint64_t)no_cells * NVERTCELL * sizeof(uint32_t) 
  	+ sizeof(uint64_t);
  fprintf( fp, "<DataArray type=\"UInt32\" Name=\"offsets\" "
  	"format=\"appended\" offset=\"%lu\"/>\n", count );
  count +=  (uint64_t)no_cells * sizeof(uint32_t) + sizeof(uint64_t);  
  fprintf( fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" "
  	"offset=\"%lu\"/>\n", count );
  count +=  (uint64_t)no_cells * sizeof(uint8_t) + sizeof(uint64_t);  
  fputs( "</Cells>\n", fp );
  fputs( "</Piece>\n", fp );
  fputs( "</UnstructuredGrid>\n", fp );
  fputs( "<AppendedData encoding=\"raw\">\n", fp );
  fputs( "_", fp );
  uint64_t block_len = (uint64_t)no_cells * sizeof(PARAVIEW_DATATYPE);
# if dimension == 2
    PARAVIEW_DATATYPE vz = 0.;
# endif
  PARAVIEW_DATATYPE fieldval = 0., vcoord = 0.;
  for (scalar s in list) 
  {
    fwrite( &block_len, sizeof(uint64_t), 1, fp );
    foreach(serial, noauto)
    {
      fieldval = (PARAVIEW_DATATYPE)(val(s));
      fwrite( &fieldval, sizeof(PARAVIEW_DATATYPE), 1, fp );
    }
  }
  block_len = no_cells * 3 * sizeof(PARAVIEW_DATATYPE);
  for (vector v in vlist) 
  {
    fwrite( &block_len, sizeof(uint64_t), 1, fp );
    foreach(serial, noauto)
    {
      fieldval = (PARAVIEW_DATATYPE)(val(v.x));
      fwrite( &fieldval, sizeof(PARAVIEW_DATATYPE), 1, fp );
      fieldval = (PARAVIEW_DATATYPE)(val(v.y));      
      fwrite( &fieldval, sizeof(PARAVIEW_DATATYPE), 1, fp );
#     if dimension == 2        
	fwrite( &vz, sizeof(PARAVIEW_DATATYPE), 1, fp );
#     endif
#     if dimension == 3
        fieldval = (PARAVIEW_DATATYPE)(val(v.z));
	fwrite( &fieldval, sizeof(PARAVIEW_DATATYPE), 1, fp );
#     endif
    }
  }
  block_len = no_points * 3 * sizeof(PARAVIEW_DATATYPE);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  foreach_vertex(serial, noauto)
  {
    vcoord = (PARAVIEW_DATATYPE)(x);
    fwrite( &vcoord, sizeof(PARAVIEW_DATATYPE), 1, fp );
    vcoord = (PARAVIEW_DATATYPE)(y);    
    fwrite( &vcoord, sizeof(PARAVIEW_DATATYPE), 1, fp );
    vcoord = (PARAVIEW_DATATYPE)(z);     
    fwrite( &vcoord, sizeof(PARAVIEW_DATATYPE), 1, fp );
  }
  // Additional duplicated vertices in case of periodicity
  PARAVIEW_DATATYPE supvertex[3];
  supvertex[2] = 0.;
  if ( Period.x || Period.y || Period.z )
  {
    foreach(serial, noauto)
      if ( per_mask[] < 0.5 )
      {
        for (size_t j=0;j<NVERTCELL;j++)
	{
	  supvertex[0] = (PARAVIEW_DATATYPE)(x + percelldir[j].x * Delta);
	  supvertex[1] = (PARAVIEW_DATATYPE)(y + percelldir[j].y * Delta);	
#         if dimension == 3
	    supvertex[2] = (PARAVIEW_DATATYPE)(z + percelldir[j].z * Delta);
#         endif	    
	  for (unsigned int k=0;k<3;k++)  
	    fwrite( &(supvertex[k]), sizeof(PARAVIEW_DATATYPE), 1, fp );
	}	  
      }
    free( percelldir );  
  }
  block_len = (uint64_t)no_cells * NVERTCELL * sizeof(uint32_t);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  uint32_t connectivity[NVERTCELL];
  foreach(serial, noauto)
  {
    if ( per_mask[] )
    {
#     if dimension == 2
        connectivity[0] = (uint32_t)(marker[]);
        connectivity[1] = (uint32_t)(marker[1,0]);
        connectivity[2] = (uint32_t)(marker[1,1]);
        connectivity[3] = (uint32_t)(marker[0,1]);
#     endif
#     if dimension == 3
        connectivity[0] = (uint32_t)(marker[]);
        connectivity[1] = (uint32_t)(marker[1,0,0]);
        connectivity[2] = (uint32_t)(marker[1,1,0]);
        connectivity[3] = (uint32_t)(marker[0,1,0]);
        connectivity[4] = (uint32_t)(marker[0,0,1]);
        connectivity[5] = (uint32_t)(marker[1,0,1]);
        connectivity[6] = (uint32_t)(marker[1,1,1]);
        connectivity[7] = (uint32_t)(marker[0,1,1]);
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        connectivity[0] = vertexnum++;
        connectivity[1] = vertexnum++;
        connectivity[2] = vertexnum++;
        connectivity[3] = vertexnum++;
#     endif
#     if dimension == 3
        connectivity[0] = vertexnum++;
        connectivity[1] = vertexnum++;
        connectivity[2] = vertexnum++;
        connectivity[3] = vertexnum++;
        connectivity[4] = vertexnum++;
        connectivity[5] = vertexnum++;
        connectivity[6] = vertexnum++;
        connectivity[7] = vertexnum++;
#     endif
    }    
    fwrite( &connectivity, sizeof(uint32_t), NVERTCELL, fp );  
  }
  block_len = (uint64_t)no_cells * sizeof(uint32_t);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  uint32_t offset = 0;
  for (uint32_t i = 1; i < no_cells+1; i++)
  {
    offset = i * NVERTCELL;
    fwrite( &offset, sizeof(uint32_t), 1, fp );
  } 
  block_len = (uint64_t)no_cells * sizeof(uint8_t);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  uint8_t ctype = CELLTYPE;
  for (uint32_t i = 1; i < no_cells+1; i++)
    fwrite( &ctype, sizeof(uint8_t), 1, fp );      
  fputs( "\n", fp );
  fputs( "</AppendedData>\n", fp );
  fputs( "</VTKFile>\n", fp );
  fflush( fp );
# if defined(_OPENMP)
    omp_set_num_threads( num_omp );
# endif
}




# if _MPI
uint8_t count_revifs( uint64_t n ) 
{
  if (n > 9999999999) return 11;
  if (n > 999999999) return 10;
  if (n > 99999999) return 9;
  if (n > 9999999) return 8;
  if (n > 999999) return 7;
  if (n > 99999) return 6;
  if (n > 9999) return 5;
  if (n > 999) return 4;
  if (n > 99) return 3;
  if (n > 9) return 2;
  return 1;
}




void output_vtu_ascii_foreach_MPIIO( scalar* list, vector* vlist )
{
  MPI_File file;
  MPI_Status status;
  MPI_Datatype num_as_string;
  char fmt[8] = "%12.5e ";
  char endfmt[8] = "%12.5e\n";
  const int charspernum = 13;
  size_t counter = 0; 
  char header[1000],line[200];
  unsigned int nprocs = npe(), i, j, m, rank = pid();  
  
  // Create a MPI datatype to write doublea as strings with a given format
  MPI_Type_contiguous( charspernum, MPI_CHAR, &num_as_string ); 
  MPI_Type_commit( &num_as_string );

  // Open the file 
  MPI_File_open( MPI_COMM_WORLD, "Res/titi.vtu", 
  	MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file );

  vertex scalar marker[];
  // We use uint64_t for the 3 following variables to avoid excessive casting 
  uint64_t no_points = 0, no_cells = 0, vertexnum = 0, length = 0,
  	vertexnumref = 0;

  // In case of periodicity
  scalar per_mask[];
  foreach(serial, noauto) per_mask[] = 1.; 
  coord* percelldir = NULL;
  if ( Period.x || Period.y || Period.z )
  {  
    percelldir = init_percelldir();
    foreach(serial, noauto)
    {
      if ( Period.x )
        if ( x + Delta > X0 + L0 || x - Delta < X0 )
          per_mask[] = 0.;
#     if dimension > 1
        if ( Period.y )
          if ( y + Delta > Y0 + L0 ||  y - Delta < Y0 )
            per_mask[] = 0.;
#     endif
#     if dimension > 2
        if ( Period.z )
          if ( z + Delta > Z0 + L0 || z - Delta < Z0 )
            per_mask[] = 0.;
#     endif
    }
  }

  foreach_vertex(serial, noauto)
  {
    marker[] = vertexnum++;
    no_points += 1;  
  }
  foreach(serial, noauto)
  {
    // Additional (duplicated) vertices per cell that have a face belonging to a
    // periodic face of the whole domain
    if ( per_mask[] < 0.5 ) no_points += NVERTCELL;
    
    // Number of cells on this process
    no_cells += 1;
  }

  
  // Points and cells on all processes
  uint64_t total_nbpts;
  MPI_Allreduce( &no_points, &total_nbpts, 1, MPI_UINT64_T, MPI_SUM, 
  	MPI_COMM_WORLD );
  uint64_t total_nbcells;
  MPI_Allreduce( &no_cells, &total_nbcells, 1, MPI_UINT64_T, MPI_SUM, 
  	MPI_COMM_WORLD );  
  uint64_t* nbpts_per_proc = (uint64_t*) calloc( nprocs, sizeof(uint64_t) );
  uint64_t* nbcells_per_proc = (uint64_t*) calloc( nprocs, sizeof(uint64_t) );
  MPI_Allgather( &no_points, 1, MPI_UINT64_T, nbpts_per_proc, 1, MPI_UINT64_T, 
  	MPI_COMM_WORLD );
  MPI_Allgather( &no_cells, 1, MPI_UINT64_T, nbcells_per_proc, 1, MPI_UINT64_T, 
  	MPI_COMM_WORLD );	


  // Header
  sprintf( header, "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
  	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n" );
  strcat( header, "<UnstructuredGrid>\n" );  
  sprintf( line, "<Piece NumberOfPoints=\"%lu\" NumberOfCells=\"%lu\">\n", 
  	total_nbpts, total_nbcells );	
  strcat( header, line );
  sprintf( line, "<Points>\n" );
  strcat( header, line );
  sprintf( line, "<DataArray type=\"%s\" "
  	"NumberOfComponents=\"3\" format=\"ascii\">\n", PARAVIEW_DATANAME );
  strcat( header, line );	
  if ( rank == 0 )
    MPI_File_write( file, header, strlen(header), MPI_CHAR, &status );  


  // Write point coordinates to the MPI file
  char* pts_coord = (char*) calloc( dimension * no_points * charspernum + 1, 
  	sizeof(char) ); 
  foreach_vertex(serial, noauto)
  {
#   if dimension == 2
      sprintf( &pts_coord[2*counter*charspernum], fmt, x );
      sprintf( &pts_coord[(2*counter+1)*charspernum], endfmt, y );      
#   endif
#   if dimension == 3
      sprintf( &pts_coord[3*counter*charspernum], fmt, x );
      sprintf( &pts_coord[(3*counter+1)*charspernum], fmt, y );
      sprintf( &pts_coord[(3*counter+2)*charspernum], endfmt, z );      
#   endif
    ++counter;
  }
  // Additional duplicated vertices in case of periodicity
  coord supvertex = {0., 0., 0.};
  if ( Period.x || Period.y || Period.z )
  {
    foreach(serial, noauto)
      if ( per_mask[] < 0.5 )
      {
#       if dimension == 2
          for (size_t j=0;j<4;j++)
	  {
	    supvertex.x = x + percelldir[j].x * Delta;
	    supvertex.y = y + percelldir[j].y * Delta;	    
            sprintf( &pts_coord[2*counter*charspernum], fmt, supvertex.x );
            sprintf( &pts_coord[(2*counter+1)*charspernum], endfmt, 
	    	supvertex.y );
	  }		
#       endif
#       if dimension == 3
          for (size_t j=0;j<8;j++)
	  {
	    supvertex.x = x + percelldir[j].x * Delta;
	    supvertex.y = y + percelldir[j].y * Delta;	
	    supvertex.z = z + percelldir[j].z * Delta;	    	    
            sprintf( &pts_coord[3*counter*charspernum], fmt, supvertex.x );
            sprintf( &pts_coord[(3*counter+1)*charspernum], fmt, supvertex.y );
            sprintf( &pts_coord[(3*counter+2)*charspernum], endfmt, 
	    	supvertex.z );
	  }	
#       endif
        ++counter;    
      }
    free( percelldir );       
  }
  
  uint64_t* mpifile_offsets = (uint64_t*) calloc( nprocs, sizeof(uint64_t) );
  mpifile_offsets[0] = strlen(header) * sizeof(char);
  for (i=1;i<nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbpts_per_proc[i-1] * dimension * charspernum * sizeof(char);

  MPI_File_write_at_all( file, mpifile_offsets[pid()], pts_coord, 
  	dimension * no_points, num_as_string, &status ); 

  free(pts_coord); pts_coord = NULL;   

    
  // Header for connectivity
  sprintf( header, "</DataArray>\n</Points>\n<Cells>\n"
  	"<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n" );
  int mpifile_offset = mpifile_offsets[nprocs-1] 
  	+ nbpts_per_proc[nprocs-1] * dimension * charspernum * sizeof(char);
  if ( rank == 0 )
    MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status ); 

  // Compute the point number shifts 
  uint64_t* shift = (uint64_t*) calloc( nprocs, sizeof(uint64_t) );
  for (m=1;m<nprocs;m++)
    for (j=0;j<m;j++)
      shift[m] += nbpts_per_proc[j];

  // Length of the connectivity string
  vertexnumref = vertexnum;
  foreach(serial, noauto)
    if ( per_mask[] )
    {
#     if dimension == 2
        uint64_t ape1 = (uint64_t)marker[] + shift[rank];
        uint64_t ape2 = (uint64_t)marker[1,0]) + shift[rank];
        uint64_t ape3 = (uint64_t)marker[1,1]) + shift[rank];
        uint64_t int ape4 = (uint64_t)marker[0,1]) + shift[rank];
	length += count_revifs( ape1 ) + count_revifs( ape2 )
		+ count_revifs( ape3 ) + count_revifs( ape4 ) + 4;		
#     endif
#     if dimension == 3
        uint64_t ape1 = (uint64_t)marker[] + shift[rank];
        uint64_t ape2 = (uint64_t)marker[1,0,0] + shift[rank];
        uint64_t ape3 = (uint64_t)marker[1,1,0] + shift[rank];
        uint64_t ape4 = (uint64_t)marker[0,1,0] + shift[rank];
        uint64_t ape5 = (uint64_t)marker[0,0,1] + shift[rank];
        uint64_t ape6 = (uint64_t)marker[1,0,1] + shift[rank];
        uint64_t ape7 = (uint64_t)marker[1,1,1] + shift[rank];
        uint64_t ape8 = (uint64_t)marker[0,1,1] + shift[rank];
	length += count_revifs( ape1 ) + count_revifs( ape2 )
		+ count_revifs( ape3 ) + count_revifs( ape4 ) 
		+ count_revifs( ape5 ) + count_revifs( ape7 )
		+ count_revifs( ape6 ) + count_revifs( ape8 ) + 8;
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        uint64_t ape1 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape2 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape3 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape4 = vertexnum + shift[rank]; vertexnum++;
	length += count_revifs( ape1 ) + count_revifs( ape2 )
		+ count_revifs( ape3 ) + count_revifs( ape4 ) + 4;    
#     endif
#     if dimension == 3
        uint64_t ape1 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape2 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape3 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape4 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape5 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape6 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape7 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape8 = vertexnum + shift[rank]; vertexnum++;
	length += count_revifs( ape1 ) + count_revifs( ape2 )
		+ count_revifs( ape3 ) + count_revifs( ape4 ) 
		+ count_revifs( ape5 ) + count_revifs( ape7 )
		+ count_revifs( ape6 ) + count_revifs( ape8 ) + 8;
#     endif
    }
  vertexnum = vertexnumref;  

  // Write connectivity
  char* connectivity = (char*) calloc( length + 1, sizeof(char) ); 
  counter = 0;  
  foreach(serial, noauto)
    if ( per_mask[] )
    {
#     if dimension == 2
        uint64_t ape1 = (uint64_t)marker[] + shift[rank];
        uint64_t ape2 = (uint64_t)marker[1,0] + shift[rank];
        uint64_t ape3 = (uint64_t)marker[1,1] + shift[rank];
        uint64_t ape4 = (uint64_t)marker[0,1] + shift[rank];
	sprintf( &connectivity[counter], "%lu ", ape1 );
        counter += count_revifs( ape1 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape2 );
        counter += count_revifs( ape2 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape3 );
        counter += count_revifs( ape3 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape4 );
        counter += count_revifs( ape4 ) + 1;		
#     endif
#     if dimension == 3
        uint64_t ape1 = (uint64_t)marker[] + shift[rank];
        uint64_t ape2 = (uint64_t)marker[1,0,0] + shift[rank];
        uint64_t ape3 = (uint64_t)marker[1,1,0] + shift[rank];
        uint64_t ape4 = (uint64_t)marker[0,1,0] + shift[rank];
        uint64_t ape5 = (uint64_t)marker[0,0,1] + shift[rank];
        uint64_t ape6 = (uint64_t)marker[1,0,1] + shift[rank];
        uint64_t ape7 = (uint64_t)marker[1,1,1] + shift[rank];
        uint64_t ape8 = (uint64_t)marker[0,1,1] + shift[rank];
	sprintf( &connectivity[counter], "%lu ", ape1 );
        counter += count_revifs( ape1 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape2 );
        counter += count_revifs( ape2 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape3 );
        counter += count_revifs( ape3 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape4 );
        counter += count_revifs( ape4 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape5 );
        counter += count_revifs( ape5 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape6 );
        counter += count_revifs( ape6 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape7 );
        counter += count_revifs( ape7 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape8 );
        counter += count_revifs( ape8 ) + 1;		
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        uint64_t ape1 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape2 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape3 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape4 = vertexnum + shift[rank]; vertexnum++;
	sprintf( &connectivity[counter], "%lu ", ape1 );
        counter += count_revifs( ape1 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape2 );
        counter += count_revifs( ape2 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape3 );
        counter += count_revifs( ape3 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape4 );
        counter += count_revifs( ape4 ) + 1;     
#     endif
#     if dimension == 3
        uint64_t ape1 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape2 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape3 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape4 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape5 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape6 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape7 = vertexnum + shift[rank]; vertexnum++;
        uint64_t ape8 = vertexnum + shift[rank]; vertexnum++;
	sprintf( &connectivity[counter], "%lu ", ape1 );
        counter += count_revifs( ape1 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape2 );
        counter += count_revifs( ape2 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape3 );
        counter += count_revifs( ape3 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape4 );
        counter += count_revifs( ape4 ) + 1; 
	sprintf( &connectivity[counter], "%lu ", ape5 );
        counter += count_revifs( ape5 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape6 );
        counter += count_revifs( ape6 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape7 );
        counter += count_revifs( ape7 ) + 1;
	sprintf( &connectivity[counter], "%lu ", ape8 );
        counter += count_revifs( ape8 ) + 1;
#     endif
    }

  uint64_t out_length = counter;
  uint64_t* out_length_per_proc = (uint64_t*) calloc( nprocs, 
  	sizeof(uint64_t) );
  MPI_Allgather( &out_length, 1, MPI_UINT64_T, out_length_per_proc, 1, 
  	MPI_UINT64_T, MPI_COMM_WORLD );
  
  mpifile_offsets[0] = mpifile_offset + strlen(header) * sizeof(char);
  for (i=1;i<nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[rank], connectivity, 
  	out_length, MPI_CHAR, &status );   

  free(connectivity); connectivity = NULL;


  // Header for offset
  sprintf( header, "\n</DataArray>\n<DataArray type=\"UInt64\" Name=\"offsets\" "
  	"format=\"ascii\">\n" );
  mpifile_offset = mpifile_offsets[nprocs-1] 
  	+ out_length_per_proc[nprocs-1] * sizeof(char);
  if ( rank == 0 )
    MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );

  // Compute the offset shifts 
  for (m=1;m<nprocs;m++)
  {
    shift[m] = 0;
    for (j=0;j<m;j++)
      shift[m] += NVERTCELL * nbcells_per_proc[j];
  }

  // Length of the offset string
  length = 0;
  for (i = 1; i < no_cells+1; i++)
    length += count_revifs( i * NVERTCELL + shift[rank] ) + 1; 

  // Write offset
  char* offset = (char*) calloc( length + 1, sizeof(char) ); 
  counter = 0;
  for (i = 1; i < no_cells+1; i++)
  {
    sprintf( &offset[counter], "%lu ", i * NVERTCELL + shift[rank] );
    counter += count_revifs( i * NVERTCELL + shift[rank] ) + 1;
  }
  
  out_length = counter;
  MPI_Allgather( &out_length, 1, MPI_INT, out_length_per_proc, 1, MPI_INT, 
  	MPI_COMM_WORLD );
  
  mpifile_offsets[0] = mpifile_offset + strlen(header) * sizeof(char);
  for (i=1;i<nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[rank], offset, 
  	out_length, MPI_CHAR, &status );   

  free(offset); offset = NULL;
    

  // Header for types
  sprintf( header, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" "
  	"format=\"ascii\">\n" );
  mpifile_offset = mpifile_offsets[nprocs-1] 
  	+ out_length_per_proc[nprocs-1] * sizeof(char);
  if ( rank == 0 )
    MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );

  // Length of the type string
  length = no_cells * dimension;

  // Write types
  char* types = (char*) calloc( length + 1, sizeof(char) ); 
  counter = 0;
  for (i = 0; i < no_cells; i++)
  {
    sprintf( &types[counter], "%hu ", CELLTYPE );
    counter += dimension;
  }
  
  out_length = counter;
  MPI_Allgather( &out_length, 1, MPI_INT, out_length_per_proc, 1, MPI_INT, 
  	MPI_COMM_WORLD );
  
  mpifile_offsets[0] = mpifile_offset + strlen(header) * sizeof(char);
  for (i=1;i<nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[rank], types, 
  	out_length, MPI_CHAR, &status );   

  free(types); types = NULL;


  // Header for field values
  sprintf( header, 
  	"\n</DataArray>\n</Cells>\n<CellData Scalars=\"scalars\">\n" );
  mpifile_offset = mpifile_offsets[nprocs-1] 
  	+ out_length_per_proc[nprocs-1] * sizeof(char);
  if ( rank == 0 )
    MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );  


  // Write for each scalar field 
  char* scalar_data = (char*) calloc( no_cells * charspernum + 1, 
  	sizeof(char) );    
  for (scalar s in list) 
  {
    mpifile_offset += strlen(header) * sizeof(char);    
    sprintf( header, "<DataArray type=\"%s\" Name=\"%s\" "
    	"format=\"ascii\">\n", PARAVIEW_DATANAME, s.name );    
    if ( rank == 0 )
      MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );  

    counter = 0;
    foreach(serial, noauto)
    {
      sprintf( &scalar_data[counter*charspernum], endfmt, val(s) );
      ++counter;
    }
    
    mpifile_offsets[0] = mpifile_offset + strlen(header) * sizeof(char);
    for (i=1;i<nprocs;i++)
      mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbcells_per_proc[i-1] * charspernum * sizeof(char);

    MPI_File_write_at_all( file, mpifile_offsets[pid()], scalar_data, 
  	no_cells, num_as_string, &status );     

    sprintf( header, "</DataArray>\n" );
    mpifile_offset = mpifile_offsets[nprocs-1] 
  	+ nbcells_per_proc[nprocs-1] * charspernum * sizeof(char);
    if ( rank == 0 )
      MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status ); 	
  }
  
  free(scalar_data); scalar_data = NULL;
    
  char* vector_data = (char*) calloc( no_cells * dimension * charspernum + 1, 
  	sizeof(char) );    
  for (vector v in vlist) 
  {
    mpifile_offset += strlen(header) * sizeof(char);    
    sprintf( header, "<DataArray type=\"%s\" "
    	"NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", 
	PARAVIEW_DATANAME, v.x.name );    
    if ( rank == 0 )
      MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );  

    counter = 0;
    foreach(serial, noauto)
    {
#     if dimension == 2
        sprintf( &vector_data[2*counter*charspernum], fmt, val(v.x) );
        sprintf( &vector_data[(2*counter+1)*charspernum], endfmt, val(v.y) );
#     endif
#     if dimension == 3
        sprintf( &vector_data[3*counter*charspernum], fmt, val(v.x) );
        sprintf( &vector_data[(3*counter+1)*charspernum], fmt, val(v.y) );
        sprintf( &vector_data[(3*counter+2)*charspernum], endfmt, val(v.z) );
#     endif
      ++counter;
    }
    
    mpifile_offsets[0] = mpifile_offset + strlen(header) * sizeof(char);
    for (i=1;i<nprocs;i++)
      mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ nbcells_per_proc[i-1] * dimension * charspernum * sizeof(char);

    MPI_File_write_at_all( file, mpifile_offsets[pid()], vector_data, 
  	dimension * no_cells, num_as_string, &status );     

    sprintf( header, "</DataArray>\n" );
    mpifile_offset = mpifile_offsets[nprocs-1] 
  	+ nbcells_per_proc[nprocs-1] * dimension * charspernum * sizeof(char);
    if ( rank == 0 )
      MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );     
  }

  free(vector_data); vector_data = NULL;
  
  // Closing text
  mpifile_offset += strlen(header) * sizeof(char);
  sprintf( header, 
  	"</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n" ); 	   
  if ( rank == 0 )
    MPI_File_write_at( file, mpifile_offset, header, strlen(header),
    	MPI_CHAR, &status );    
  
  
  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  MPI_Type_free( &num_as_string ); 
  free(nbpts_per_proc); nbpts_per_proc = NULL;
  free(nbcells_per_proc); nbcells_per_proc = NULL;
  free(mpifile_offsets); mpifile_offsets = NULL; 
  free(out_length_per_proc); out_length_per_proc = NULL;    
}




// ----------------------------------------------------------------------------
// Methods to write binary data
void check_allocated_binary( int size )  
{
  if ( OFFSET + size >= ALLOCATED ) 
  {
    int new_size = max( 2*ALLOCATED, (int)1024 ) ;
    new_size = max( new_size, 2*(OFFSET+size) ) ;
    new_size = 4 * ( new_size/4 +1 ) ; // alignment on 4 bytes
      
    char* new_buffer = (char*) calloc( new_size, sizeof(char) ); 
    for( int i=0 ;i<OFFSET ;i++ ) new_buffer[i] = BUFFER[i] ;
    if ( BUFFER != 0 ) free(BUFFER) ;
    BUFFER = new_buffer ;
    ALLOCATED = new_size ;      
  }
}




int store_int_binary( int val )  
{
  int result = OFFSET ;
  *((int*)&(BUFFER[OFFSET])) = val ;
  OFFSET += sizeof_Int32  ;
  return result ;
}




void start_output_binary( int size, int number )
{
  int current_output_size = size*number ;
  int ncomp = current_output_size + (current_output_size+999)/1000 
  	+ 12 + sizeof_Int32 ;	
  check_allocated_binary( ncomp ) ;
  CURRENT_LENGTH = store_int_binary( current_output_size ) ;
}




void write_double_binary( double val )  
{
  *((float*)&(BUFFER[OFFSET])) = (float)val ;
  OFFSET += sizeof_Float32  ;
}




void write_int_binary( int val )  
{
//  store_int_binary(val) ;
  *((int*)&(BUFFER[OFFSET])) = val ;
  OFFSET += sizeof_Int32  ;
}




void compress_segment_binary( int seg )  
{
   static int BlockSize = 32768 ;
   int size = (int)(*((int*)&BUFFER[seg])), i ;
   
   int numFullBlocks = size / BlockSize;
   int lastBlockSize = size % BlockSize;
   int numBlocks = numFullBlocks + (lastBlockSize?1:0);

   int headerLength = numBlocks+3;

   int* CompressionHeader = (int*) calloc( headerLength, sizeof(int) ); 
   CompressionHeader[0] = numBlocks;
   CompressionHeader[1] = BlockSize;
   CompressionHeader[2] = lastBlockSize;

   unsigned long encoded_buff_size = max(BlockSize,size)  ;
   unsigned char* encoded_buff = (unsigned char*) calloc( encoded_buff_size, 
   	sizeof(unsigned char) ); 
   int encoded_offset = 0, block ;
   for( block=0 ; block<numBlocks ; block++ )
   {
      int buffer_start = seg + sizeof_Int32 + block*BlockSize ;
      int length = ( block+1<numBlocks || !lastBlockSize ? 
      	BlockSize : lastBlockSize ) ;
      unsigned char* to_encode = (unsigned char *)(&BUFFER[buffer_start]) ;
      unsigned char* encoded = &encoded_buff[encoded_offset] ;
      unsigned long ncomp = encoded_buff_size - encoded_offset ;

      compress2( (Bytef*)encoded, &ncomp, (const Bytef*)to_encode,
	length, Z_DEFAULT_COMPRESSION );
	 
      CompressionHeader[3+block] = (int)(ncomp) ;
      encoded_offset += (int)(ncomp) ;      
   }
   
   OFFSET = seg ;
   check_allocated_binary( headerLength * sizeof_Int32 + encoded_offset ) ;
   
   for( i=0 ; i<headerLength ; i++ )
      store_int_binary( CompressionHeader[i] ) ;     

   for( i=0 ; i<encoded_offset ; i++ )
      BUFFER[OFFSET++] = encoded_buff[i] ;

   if( OFFSET%4 != 0 )
      OFFSET = 4*( OFFSET/4 +1 ) ; // Re-alignment
   
   free(CompressionHeader); CompressionHeader = NULL;
   free(encoded_buff); encoded_buff = NULL;
}




void flush_binary()  
{
  compress_segment_binary( CURRENT_LENGTH ) ;         
}
// End of Methods to write binary data
// ----------------------------------------------------------------------------




void output_vtu_bin_foreach_MPIIO( scalar* list, vector* vlist )
{
  MPI_File file;
  MPI_Status status;
  int point_binary_offset = 0, connectivity_binary_offset, 
  	offsets_binary_offset, cellstype_binary_offset, total_offset;
  char header[5000], line[500];
  int nprocs = npe(), i, m, rank = pid(); 
  	  
  // Open the file 
  MPI_File_open( MPI_COMM_WORLD, "Res/titi_bin.vtu", 
  	MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file );  


  vertex scalar marker[];
  // We use unsigned int for the 2 following variables because the number
  // of cells/vertices on a single proc never exceeds the limit of 4,294,967,295
  // on 64-bit systems. Hence there is not need to use long unsigned int
  unsigned int no_points = 0, no_cells = 0;
  int vertexnum = 0;

  // In case of periodicity
  scalar per_mask[];
  foreach(serial, noauto) per_mask[] = 1.; 
  coord* percelldir = NULL;
  if ( Period.x || Period.y || Period.z )
  {  
    percelldir = init_percelldir();
    foreach(serial, noauto)
    {
      if ( Period.x )
        if ( x + Delta > X0 + L0 || x - Delta < X0 )
          per_mask[] = 0.;
#     if dimension > 1
        if ( Period.y )
          if ( y + Delta > Y0 + L0 ||  y - Delta < Y0 )
            per_mask[] = 0.;
#     endif
#     if dimension > 2
        if ( Period.z )
          if ( z + Delta > Z0 + L0 || z - Delta < Z0 )
            per_mask[] = 0.;
#     endif
    }
  }

  foreach_vertex(serial, noauto)
  {
    marker[] = vertexnum++;
    no_points += 1;  
  }
  foreach(serial, noauto)
  {
    // Additional (duplicated) vertices per cell that have a face belonging to a
    // periodic face of the whole domain
    if ( per_mask[] < 0.5 ) no_points += NVERTCELL;
    
    // Number of cells on this process
    no_cells += 1;
  }


  // Write point coordinates to the binary buffer
  start_output_binary( sizeof_Float32, dimension * no_points ) ;
  foreach_vertex(serial, noauto)
  {
#    if dimension == 2
       write_double_binary( x ) ;
       write_double_binary( y ) ;     
#    endif
#    if dimension == 3
       write_double_binary( x ) ;
       write_double_binary( y ) ;     
       write_double_binary( z ) ;
#    endif       
  }     
  // Additional duplicated vertices in case of periodicity
  coord supvertex = {0., 0., 0.};
  if ( Period.x || Period.y || Period.z )
  {
    foreach(serial, noauto)
      if ( per_mask[] < 0.5 )
      {
#       if dimension == 2
          for (size_t j=0;j<4;j++)
	  {
	    supvertex.x = x + percelldir[j].x * Delta;
	    supvertex.y = y + percelldir[j].y * Delta;
	    write_double_binary( supvertex.x ) ;
	    write_double_binary( supvertex.y ) ;	    
	  }		
#       endif
#       if dimension == 3
          for (size_t j=0;j<8;j++)
	  {
	    supvertex.x = x + percelldir[j].x * Delta;
	    supvertex.y = y + percelldir[j].y * Delta;	
	    supvertex.z = z + percelldir[j].z * Delta;	    	    
	    write_double_binary( supvertex.x ) ;
	    write_double_binary( supvertex.y ) ;
	    write_double_binary( supvertex.z ) ;	    
	  }	
#       endif
      }
    free( percelldir );       
  }
  compress_segment_binary( CURRENT_LENGTH );    


  // Write connectivity to the binary buffer
  connectivity_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, NVERTCELL * no_cells ) ;  
  foreach(serial, noauto)
    if ( per_mask[] )
    {
#     if dimension == 2
        int ape1 = (int)marker[];
        int ape2 = (int)marker[1,0];
        int ape3 = (int)marker[1,1];
        int ape4 = (int)marker[0,1];
	write_int_binary( ape1 );
        write_int_binary( ape2 );
	write_int_binary( ape3 );
        write_int_binary( ape4 );				
#     endif
#     if dimension == 3
        int ape1 = (int)marker[];
        int ape2 = (int)marker[1,0,0];
        int ape3 = (int)marker[1,1,0];
        int ape4 = (int)marker[0,1,0];
        int ape5 = (int)marker[0,0,1];
        int ape6 = (int)marker[1,0,1];
        int ape7 = (int)marker[1,1,1];
        int ape8 = (int)marker[0,1,1];
	write_int_binary( ape1 );
        write_int_binary( ape2 );
	write_int_binary( ape3 );
        write_int_binary( ape4 );
	write_int_binary( ape5 );
        write_int_binary( ape6 );
	write_int_binary( ape7 );
        write_int_binary( ape8 );			
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        int ape1 = vertexnum; vertexnum++;
        int ape2 = vertexnum; vertexnum++;
        int ape3 = vertexnum; vertexnum++;
        int ape4 = vertexnum; vertexnum++;
	write_int_binary( ape1 );
        write_int_binary( ape2 );
	write_int_binary( ape3 );
        write_int_binary( ape4 );   
#     endif
#     if dimension == 3
        int ape1 = vertexnum; vertexnum++;
        int ape2 = vertexnum; vertexnum++;
        int ape3 = vertexnum; vertexnum++;
        int ape4 = vertexnum; vertexnum++;
        int ape5 = vertexnum; vertexnum++;
        int ape6 = vertexnum; vertexnum++;
        int ape7 = vertexnum; vertexnum++;
        int ape8 = vertexnum; vertexnum++;
	write_int_binary( ape1 );
        write_int_binary( ape2 );
	write_int_binary( ape3 );
        write_int_binary( ape4 );
	write_int_binary( ape5 );
        write_int_binary( ape6 );
	write_int_binary( ape7 );
        write_int_binary( ape8 );
#     endif
    }  
  compress_segment_binary( CURRENT_LENGTH );  


  // Write offsets to the binary buffer  
  offsets_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, no_cells ) ;  
  for (i = 1; i < no_cells+1; i++) 
    write_int_binary( i * NVERTCELL );
  compress_segment_binary( CURRENT_LENGTH ); 
  
  
  // Write cell types to the binary buffer  
  cellstype_binary_offset = OFFSET;
  start_output_binary( sizeof_Int32, no_cells ) ;  
  for (i = 0; i < no_cells; i++)
    write_int_binary( CELLTYPE );  
  compress_segment_binary( CURRENT_LENGTH );   


  // Scalar field values
  int nscalar = 0, iscal = 0;
  for (scalar s in list) ++nscalar;
  int* scalar_binary_offset = (int*) calloc( nscalar, sizeof(int) );

  // Write each scalar field to the binary buffer 
  for (scalar s in list) 
  {
    scalar_binary_offset[iscal] = OFFSET;
    start_output_binary( sizeof_Float32, no_cells );    
    foreach(serial, noauto) write_double_binary( val(s) );
    compress_segment_binary( CURRENT_LENGTH );
    ++iscal;      
  }
  
  
  // Vector field values
  int nvector = 0, ivec = 0;
  for (vector v in vlist)  ++nvector;
  int* vector_binary_offset = (int*) calloc( nvector, sizeof(int) );

  // Write each scalar field to the binary buffer 
  for (vector v in vlist)
  {
    vector_binary_offset[ivec] = OFFSET;
    start_output_binary( sizeof_Float32, dimension * no_cells );     
    foreach(serial, noauto)
    {
#     if dimension == 2
        write_double_binary( val(v.x) );
        write_double_binary( val(v.y) );
#     endif
#     if dimension == 3
        write_double_binary( val(v.x) );
        write_double_binary( val(v.y) );
        write_double_binary( val(v.z) );      
#     endif
    }    
    compress_segment_binary( CURRENT_LENGTH );
    ++ivec;      
  }    


  total_offset = OFFSET; 
  int* total_binary_offset_per_proc = (int*) calloc( nprocs, sizeof(int) );
  MPI_Allgather( &total_offset, 1, MPI_INT, total_binary_offset_per_proc, 1, 
  	MPI_INT, MPI_COMM_WORLD );    
  int* cumul_binary_offset_per_proc = (int*) calloc( nprocs, sizeof(int) );
  cumul_binary_offset_per_proc[0] = 0;
  for (m=1;m<nprocs;m++) 
    cumul_binary_offset_per_proc[m] = cumul_binary_offset_per_proc[m-1]
    	+ total_binary_offset_per_proc[m-1];


  // Header per piece + general for 1st and last process
  if ( rank == 0 )
  {   
    sprintf( header, "<?xml version=\"1.0\"?>\n"
    	"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	"byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n"
    	"<UnstructuredGrid>\n<Piece NumberOfPoints=\"%u\" "
	"NumberOfCells=\"%u\">\n", no_points, no_cells );  
  }
  else
  {
    sprintf( header, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n", 
  	no_points, no_cells );
  }	
  sprintf( line, "<Points>\n" );
  strcat( header, line );  
  sprintf( line, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
  	"offset=\"%d\" format=\"appended\"></DataArray>\n</Points>\n", 
	cumul_binary_offset_per_proc[rank] + point_binary_offset );
  strcat( header, line );     
  sprintf( line, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\"" 
  	" offset=\"%d\" format=\"appended\"></DataArray>\n",	
	cumul_binary_offset_per_proc[rank] + connectivity_binary_offset );
  strcat( header, line );  	
  sprintf( line, "<DataArray type=\"Int32\" Name=\"offsets\" offset=\"%d"
  	"\" format=\"appended\"></DataArray>\n",	
	cumul_binary_offset_per_proc[rank] + offsets_binary_offset );
  strcat( header, line ); 		
  sprintf( line, "<DataArray type=\"Int32\" Name=\"types\" offset=\"%d"
  	"\" format=\"appended\"></DataArray>\n</Cells>\n",
	cumul_binary_offset_per_proc[rank] + cellstype_binary_offset );   
  strcat( header, line );
  sprintf( line, "<CellData Scalars=\"scalars\">\n" );
  strcat( header, line );  
  iscal = 0;
  for (scalar s in list) 
  {
    sprintf( line, "<DataArray type=\"Float32\" Name=\"%s\" offset=\""
    	"%d\" format=\"appended\"></DataArray>\n", 
	s.name, cumul_binary_offset_per_proc[rank] 
	+ scalar_binary_offset[iscal] );
    strcat( header, line );  	
    ++iscal;  
  }
  ivec = 0;
  for (vector v in vlist)
  {  
    sprintf( line, "<DataArray type=\"Float32\" Name=\"%s\" "
    	"NumberOfComponents=\"3\" offset=\""
    	"%d\" format=\"appended\"></DataArray>\n", 
	v.x.name, cumul_binary_offset_per_proc[rank] 
	+ vector_binary_offset[ivec] );
    strcat( header, line );  	  
    ++ivec;
  }
  sprintf( line, "</CellData>\n</Piece>\n" );
  strcat( header, line );   
  if ( rank == nprocs - 1 )
  { 
    sprintf( line, 
    	"</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_" );
    strcat( header, line ); 	  
  }            


  int header_length = strlen( header );
  int* header_length_per_proc = (int*) calloc( nprocs, sizeof(int) );  
  MPI_Allgather( &header_length, 1, MPI_INT, header_length_per_proc, 1, 
  	MPI_INT, MPI_COMM_WORLD ); 

  int* mpifile_offsets = (int*) calloc( nprocs, sizeof(int) );
  mpifile_offsets[0] = 0;
  for (m=1;m<nprocs;m++) 
    mpifile_offsets[m] = mpifile_offsets[m-1]
    	+ header_length_per_proc[m-1] * sizeof(char);  
  MPI_File_write_at_all( file, mpifile_offsets[rank], header, 
  	header_length, MPI_CHAR, &status ); 


  // Write the binary buffers
  int starting_point = mpifile_offsets[nprocs-1] 
  	+ header_length_per_proc[nprocs-1] * sizeof(char); 
  for (m=0;m<nprocs;m++) 
    mpifile_offsets[m] = starting_point
    	+ cumul_binary_offset_per_proc[m] * sizeof(char);	
  MPI_File_write_at_all( file, mpifile_offsets[rank], BUFFER, OFFSET, 
	MPI_CHAR, &status );	 	
  free(BUFFER) ; BUFFER = NULL ;
  ALLOCATED = 0 ;
  OFFSET = 0 ;


  // Closing text
  sprintf( header, "\n</AppendedData>\n</VTKFile>\n" );
  header_length = strlen( header );
  starting_point = mpifile_offsets[nprocs-1] 
  	+ total_binary_offset_per_proc[nprocs-1] * sizeof(char);
  if ( rank == 0 )
    MPI_File_write_at( file, starting_point, header, header_length, 
    	MPI_CHAR, &status );

  
  // Close MPI file and free remaining pointers
  MPI_File_close( &file );
  free(total_binary_offset_per_proc); total_binary_offset_per_proc = NULL;
  free(cumul_binary_offset_per_proc); cumul_binary_offset_per_proc = NULL;
  free(header_length_per_proc); header_length_per_proc = NULL;
  free(mpifile_offsets); mpifile_offsets = NULL; 
  free(scalar_binary_offset); scalar_binary_offset = NULL;      
}
# endif
