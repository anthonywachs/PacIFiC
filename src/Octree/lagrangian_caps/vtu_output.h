/**
# Wrapper for output functions with Paraview

File copy-pasted (with minor modifications by Can Celçuk and Guodong Gai) from [Oystein Lande's sandbox](http://basilisk.fr/sandbox/oystelan/output_vtu_foreach.h). All credit goes to Oystein Lande, Can Selçuk and Guodong Gai.

*/

/** 
# Wrapper for output functions with Paraview
*/

# ifndef PARAVIEW_DATATYPE_DOUBLE
#   define PARAVIEW_DATATYPE_DOUBLE 0 // 1 for double and 0 for float
# endif
# if ( PARAVIEW_DATATYPE_DOUBLE == 1 )
#   define PARAVIEW_DATATYPE double
#   define PARAVIEW_DATANAME "Float64"
# else
#   define PARAVIEW_DATATYPE float
#   define PARAVIEW_DATANAME "Float32"
# endif
# ifndef PARAVIEW_BINFILE
#   define PARAVIEW_BINFILE 1
# endif

// # include "DLMFD_Output_vtu_foreach.h"


/*Here we defined the directory of the results*/
# define result_dir "Res" 
# define result_fluid_rootfilename "fluid_basilisk"



/** 
# Paraview functions 
*/

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

  for (int i = 0; i < npe(); i++)
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
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
Bug correction: %g turns into scientific notation for high integer values. This 
is not supported by paraview. Hence a fix was needed. Oystein Lande 2017
*/
void output_vtu_ascii_foreach( scalar* list, vector* vlist, 
	FILE* fp, bool linear )
{

# if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads( 1 );
# endif

  vertex scalar marker[];
  // We use unsigned int for the 3 following variables because the number
  // of cells/vertices on a single proc never exceeds the limit of 4,294,967,295
  // on 64-bit systems. Hence there is not need to use long unsigned int
  unsigned int no_points = 0, no_cells = 0, vextexnum = 0 ;

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
        if ( x + Delta > X0 + L0*L0_ratio.x || x - Delta < X0 )
          per_mask[] = 0.;
#     if dimension > 1
        if ( Period.y )
          if ( y + Delta > Y0 + L0*L0_ratio.y ||  y - Delta < Y0 )
            per_mask[] = 0.;
#     endif
#     if dimension > 2
        if ( Period.z )
          if ( z + Delta > Z0 + L0*L0_ratio.z || z - Delta < Z0 )
            per_mask[] = 0.;
#     endif
    }
  }

  foreach_vertex(serial, noauto)
  {
    marker[] = vextexnum++;
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
  fputs( "<DataArray type=\"Int64\" Name=\"connectivity\" "
  	"format=\"ascii\">\n", fp );
  foreach(serial, noauto)
    if ( per_mask[] )
    {
#     if dimension == 2
        unsigned int ape1 = marker[];
        unsigned int ape2 = marker[1,0];
        unsigned int ape3 = marker[1,1];
        unsigned int ape4 = marker[0,1];
        fprintf( fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4 );
#     endif
#     if dimension == 3
        unsigned int ape1 = marker[];
        unsigned int ape2 = marker[1,0,0];
        unsigned int ape3 = marker[1,1,0];
        unsigned int ape4 = marker[0,1,0];
        unsigned int ape5 = marker[0,0,1];
        unsigned int ape6 = marker[1,0,1];
        unsigned int ape7 = marker[1,1,1];
        unsigned int ape8 = marker[0,1,1];
        fprintf( fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, 
      		ape6, ape7, ape8 );
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        unsigned int ape1 = vextexnum++;
        unsigned int ape2 = vextexnum++;
        unsigned int ape3 = vextexnum++;
        unsigned int ape4 = vextexnum++;
        fprintf( fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4 );
#     endif
#     if dimension == 3
        unsigned int ape1 = vextexnum++;
        unsigned int ape2 = vextexnum++;
        unsigned int ape3 = vextexnum++;
        unsigned int ape4 = vextexnum++;
        unsigned int ape5 = vextexnum++;
        unsigned int ape6 = vextexnum++;
        unsigned int ape7 = vextexnum++;
        unsigned int ape8 = vextexnum++;
        fprintf( fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, 
      		ape6, ape7, ape8 );
#     endif
    }
  fputs( "</DataArray>\n", fp );
  fputs( "<DataArray type=\"Int64\" Name=\"offsets\" "
  	"format=\"ascii\">\n", fp );
  for (unsigned int i = 1; i < no_cells+1; i++)
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
output_pvtu_bin() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
Bug correction: %g turns into scientific notation for high integer values. 
This is not supported by paraview. Hence a fix was needed.
Oystein Lande 2017
*/
void output_vtu_bin_foreach( scalar* list, vector* vlist, FILE* fp, 
	bool linear )
{

# if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads( 1 );
# endif

  vertex scalar marker[];
  // We use unsigned int for the 3 following variables because the number
  // of cells/vertices on a single proc never exceeds the limit of 4,294,967,295
  // on 64-bit systems. Hence there is not need to use long unsigned int  
  unsigned int no_points = 0, no_cells = 0, vextexnum = 0 ;

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
        if ( x + Delta > X0 + L0*L0_ratio.x || x - Delta < X0 )
          per_mask[] = 0.;
#     if dimension > 1
        if ( Period.y )
          if ( y + Delta > Y0 + L0*L0_ratio.y ||  y - Delta < Y0 )
            per_mask[] = 0.;
#     endif
#     if dimension > 2
        if ( Period.z )
          if ( z + Delta > Z0 + L0*L0_ratio.z || z - Delta < Z0 )
            per_mask[] = 0.;
#     endif
    }
  }

  foreach_vertex(serial, noauto)
  {
    marker[] = vextexnum++;
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
    count += no_cells * sizeof(PARAVIEW_DATATYPE) + sizeof(uint64_t);    
    fputs( "</DataArray>\n", fp );
  }
  for (vector v in vlist) 
  {
    fprintf( fp,"<DataArray type=\"%s\" Name=\"%s\" "
    	"NumberOfComponents=\"3\" format=\"appended\" offset=\"%lu\">\n", 
	PARAVIEW_DATANAME, v.x.name, count );
    count += 3 * no_cells * sizeof(PARAVIEW_DATATYPE) + sizeof(uint64_t);
    fputs( "</DataArray>\n", fp );
  }
  fputs( "</CellData>\n", fp );
  fputs( "<Points>\n", fp );
  fprintf( fp, "<DataArray type=\"%s\" NumberOfComponents=\"3\" "
  	"format=\"appended\" offset=\"%lu\">\n", PARAVIEW_DATANAME, count );
  count += 3 * no_points * sizeof(PARAVIEW_DATATYPE) + sizeof(uint64_t);
  fputs( "</DataArray>\n", fp );
  fputs( "</Points>\n", fp );
  fputs( "<Cells>\n", fp );
  fprintf( fp,"<DataArray type=\"UInt32\" Name=\"connectivity\" "
  	"format=\"appended\" offset=\"%lu\"/>\n", count );
  count +=  no_cells * NVERTCELL * sizeof(unsigned int) + sizeof(uint64_t);
  fprintf( fp, "<DataArray type=\"UInt32\" Name=\"offsets\" "
  	"format=\"appended\" offset=\"%lu\"/>\n", count );
  count +=  no_cells * sizeof(unsigned int) + sizeof(uint64_t);  
  fprintf( fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" "
  	"offset=\"%lu\"/>\n", count );
  count +=  no_cells * sizeof(uint8_t) + sizeof(uint64_t);  
  fputs( "</Cells>\n", fp );
  fputs( "</Piece>\n", fp );
  fputs( "</UnstructuredGrid>\n", fp );
  fputs( "<AppendedData encoding=\"raw\">\n", fp );
  fputs( "_", fp );
  uint64_t block_len = no_cells * sizeof(PARAVIEW_DATATYPE);
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
	  supvertex[2] = (PARAVIEW_DATATYPE)(z + percelldir[j].z * Delta);
	  for (unsigned int k=0;k<3;k++)
	    fwrite( &(supvertex[k]), sizeof(PARAVIEW_DATATYPE), 1, fp );
	}	  
      }
    free( percelldir );  
  }
  block_len = no_cells * NVERTCELL * sizeof(unsigned int);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  unsigned int connectivity[NVERTCELL];
  foreach(serial, noauto)
  {
    if ( per_mask[] )
    {
#     if dimension == 2
        connectivity[0] = (unsigned int)(marker[]);
        connectivity[1] = (unsigned int)(marker[1,0]);
        connectivity[2] = (unsigned int)(marker[1,1]);
        connectivity[3] = (unsigned int)(marker[0,1]);
#     endif
#     if dimension == 3
        connectivity[0] = (unsigned int)(marker[]);
        connectivity[1] = (unsigned int)(marker[1,0,0]);
        connectivity[2] = (unsigned int)(marker[1,1,0]);
        connectivity[3] = (unsigned int)(marker[0,1,0]);
        connectivity[4] = (unsigned int)(marker[0,0,1]);
        connectivity[5] = (unsigned int)(marker[1,0,1]);
        connectivity[6] = (unsigned int)(marker[1,1,1]);
        connectivity[7] = (unsigned int)(marker[0,1,1]);
#     endif
    }
    // Additional duplicated vertices
    else
    {
#     if dimension == 2
        connectivity[0] = vextexnum++;
        connectivity[1] = vextexnum++;
        connectivity[2] = vextexnum++;
        connectivity[3] = vextexnum++;
#     endif
#     if dimension == 3
        connectivity[0] = vextexnum++;
        connectivity[1] = vextexnum++;
        connectivity[2] = vextexnum++;
        connectivity[3] = vextexnum++;
        connectivity[4] = vextexnum++;
        connectivity[5] = vextexnum++;
        connectivity[6] = vextexnum++;
        connectivity[7] = vextexnum++;
#     endif
    }    
    fwrite( &connectivity, sizeof(unsigned int), NVERTCELL, fp );  
  }
  block_len = no_cells * sizeof(unsigned int);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  unsigned int offset = 0;
  for (unsigned int i = 1; i < no_cells+1; i++)
  {
    offset = i * NVERTCELL;
    fwrite( &offset, sizeof(unsigned int), 1, fp );
  } 
  block_len = no_cells * sizeof(uint8_t);
  fwrite( &block_len, sizeof(uint64_t), 1, fp );
  int8_t ctype = CELLTYPE;
  for (unsigned int i = 1; i < no_cells+1; i++)
    fwrite( &ctype, sizeof(int8_t), 1, fp );      
  fputs( "\n", fp );
  fputs( "</AppendedData>\n", fp );
  fputs( "</VTKFile>\n", fp );
  fflush( fp );
# if defined(_OPENMP)
    omp_set_num_threads( num_omp );
# endif
}

//----------------------------------------------------------------------------
void output_pvd( FILE * fp, char const* times_series )
//----------------------------------------------------------------------------
{
  fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"Collection\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp );
  fputs( "<Collection>\n", fp );
  fputs( times_series, fp );
  fputs( "</Collection>\n", fp );
  fputs( "</VTKFile>\n", fp );
}



//----------------------------------------------------------------------------
void save_data( scalar* list, vector* vlist, double const time , int cycle_number)
//----------------------------------------------------------------------------
{
  char vtk_times_series[100000] = "";
  // static int cycle_number = 0;
  //if ( !cycle_number ) cycle_number = init_cycle_number;//ggd
  
  FILE * fpvtk;
  char filename_vtu[80] = "";
  char filename_pvtu[80] = "";     
  char suffix[80] = "";

  // Write the VTU file
  sprintf( filename_vtu, "%s", result_dir );
  strcat( filename_vtu, "/" );  
  strcat( filename_vtu, result_fluid_rootfilename );
  
  sprintf( suffix, "_T%d_%d.vtu", cycle_number, pid() );
  strcat( filename_vtu, suffix );
 
  fpvtk = fopen( filename_vtu, "w" );
  if ( PARAVIEW_BINFILE ) output_vtu_bin_foreach( list, vlist, fpvtk, false );
  else output_vtu_ascii_foreach( list, vlist, fpvtk, false );  
  fclose( fpvtk );
   
  // Write the PVTU file  
  if ( pid() == 0 ) 
  {
    sprintf( filename_pvtu, "%s", result_dir );
    strcat( filename_pvtu, "/" );  
    strcat( filename_pvtu, result_fluid_rootfilename );    
    sprintf( suffix, "_T%d.pvtu", cycle_number );
    strcat( filename_pvtu, suffix );

    fpvtk = fopen( filename_pvtu, "w" );
    
    sprintf( filename_vtu, "%s", result_fluid_rootfilename );
    sprintf( suffix, "_T%d", cycle_number );
    strcat( filename_vtu, suffix );
    if ( PARAVIEW_BINFILE ) output_pvtu_bin( list, vlist, fpvtk, filename_vtu );
    else output_pvtu_ascii( list, vlist, fpvtk, filename_vtu );    

    fclose( fpvtk );
  }
  
  // Write the PVD file  
  if ( pid() == 0 ) 
  {  
    char filename_pvd[80] = "";
    sprintf( filename_pvd, "%s", result_dir );
    strcat( filename_pvd, "/" );  
    strcat( filename_pvd, result_fluid_rootfilename );
    strcat( filename_pvd, ".pvd" ); 

    fpvtk = fopen( filename_pvd, "w" );

    char time_line[200] = "";
    strcpy( time_line, "<DataSet timestep=" );
    sprintf( suffix, "\"%.4e\"", time );
    strcat( time_line, suffix );
    strcat( time_line, " group=\"\" part=\"0\" file=\"" );
    strcpy( filename_pvtu, result_fluid_rootfilename );    
    sprintf( suffix, "_T%d.pvtu", cycle_number );
    strcat( filename_pvtu, suffix );
    strcat( time_line, filename_pvtu );        
    strcat( time_line, "\"/>\n" );  
    strcat( vtk_times_series, time_line );    
    output_pvd( fpvtk, vtk_times_series );
  
    fclose( fpvtk );
  }
  
  // Write the last cycle number in a file for restart  
  if ( pid() == 0 ) 
  {
    char filename_lcn[256] = ""; 
    sprintf( filename_lcn, "%s", result_dir );
    strcat( filename_lcn, "/" );  
    strcat( filename_lcn, result_fluid_rootfilename );
    strcat( filename_lcn, "_lcn_vtk.txt" );

    fpvtk = fopen( filename_lcn, "w" );    

    fprintf( fpvtk, "%d\n", cycle_number );
    
    fclose( fpvtk );            
  }
  
  ++cycle_number;        
}



//----------------------------------------------------------------------------
void save_data_ggd( scalar * list, vector * vlist, double const time , int cycle_number)
//----------------------------------------------------------------------------
{

  char vtk_times_series[100000] = "";
  // static int cycle_number = 0;
  //if ( !cycle_number ) cycle_number = init_cycle_number;//ggd

  FILE * fpvtk;
  char filename_vtu[80] = "";
  char filename_pvtu[80] = "";
  char suffix[80] = "";

  // Write the VTU file
  sprintf( filename_vtu, "%s", result_dir );
  strcat( filename_vtu, "/" );
  strcat( filename_vtu, result_fluid_rootfilename );

  sprintf( suffix, "_T%d_%d.vtu", cycle_number, pid());
  strcat( filename_vtu, suffix );

  fpvtk = fopen( filename_vtu, "w" );
  output_vtu_bin_foreach( list, vlist, fpvtk, false );
  fclose( fpvtk );

  // Write the PVTU file
  if ( pid() == 0 )
  {
    sprintf( filename_pvtu, "%s", result_dir );
    strcat( filename_pvtu, "/" );
    strcat( filename_pvtu, result_fluid_rootfilename );
    sprintf( suffix, "_T%d.pvtu", cycle_number );
    strcat( filename_pvtu, suffix );

    fpvtk = fopen( filename_pvtu, "w" );

    sprintf ( filename_vtu, "%s", result_fluid_rootfilename );
    sprintf( suffix, "_T%d", cycle_number );
    strcat( filename_vtu, suffix );
    output_pvtu_bin ( list, vlist, fpvtk, filename_vtu );

    fclose( fpvtk );
  }

  // Write the PVD file
  if ( pid() == 0 )
  {
    char filename_pvd[80] = "";
    sprintf( filename_pvd, "%s", result_dir );
    strcat( filename_pvd, "/" );
    strcat( filename_pvd, result_fluid_rootfilename );
    strcat( filename_pvd, ".pvd" );

    fpvtk = fopen( filename_pvd, "w" );

    char time_line[200] = "";
    strcpy( time_line, "<DataSet timestep=" );
    sprintf( suffix, "\"%.4e\"", time );
    strcat( time_line, suffix );
    strcat( time_line, " group=\"\" part=\"0\" file=\"" );
    strcpy( filename_pvtu, result_fluid_rootfilename );
    sprintf( suffix, "_T%d.pvtu", cycle_number );
    strcat( filename_pvtu, suffix );
    strcat( time_line, filename_pvtu );
    strcat( time_line, "\"/>\n" );
    strcat( vtk_times_series, time_line );
    output_pvd( fpvtk, vtk_times_series );

    fclose( fpvtk );
  }

  // Write the last cycle number in a file for restart
  if ( pid() == 0 )
  {
    char filename_lcn[256] = "";
    sprintf( filename_lcn, "%s", result_dir );
    strcat( filename_lcn, "/" );
    strcat( filename_lcn, result_fluid_rootfilename );
    strcat( filename_lcn, "_lcn_vtk.txt" );

    fpvtk = fopen( filename_lcn, "w" );

    fprintf( fpvtk, "%d\n", cycle_number );

    fclose( fpvtk );
  }

  ++cycle_number;
}



// //----------------------------------------------------------------------------
// void reinitialize_vtk_restart( void )
// //----------------------------------------------------------------------------
// {
//   // Get the last cycle cumber from previous simulation
//   char filename_lcn[80] = "";
//   sprintf( filename_lcn, "%s", result_dir );
//   strcat( filename_lcn, "/" );
//   strcat( filename_lcn, result_fluid_rootfilename );
//   strcat( filename_lcn, "_lcn_vtk.txt" );

//   FILE * fpvtk = fopen( filename_lcn, "r" );

//   fscanf ( fpvtk, "%d", &init_cycle_number );
//   ++init_cycle_number;

//   fclose( fpvtk );
// }


/*
//----------------------------------------------------------------------------
void reinitialize_vtk_restart( void )
//----------------------------------------------------------------------------
{
  // Get the last cycle cumber from previous simulation
  char filename_lcn[80] = "";
  sprintf( filename_lcn, "%s", result_dir );
  strcat( filename_lcn, "/" );
  strcat( filename_lcn, result_fluid_rootfilename );
  strcat( filename_lcn, "_lcn_vtk.txt" );

  FILE * fpvtk = fopen( filename_lcn, "r" );

  fscanf ( fpvtk, "%d", &init_cycle_number );
  ++init_cycle_number;

  fclose( fpvtk );

  // Re-initialize the time output series string
  if ( pid() == 0 )
  {
    char filename_pvd[80] = "";
    char time_line[256] = "";
    char start[20] = "";
    char start_ref[20] = "<DataSet";
    sprintf( filename_pvd, "%s", result_dir );
    strcat( filename_pvd, "/" );
    strcat( filename_pvd, result_fluid_rootfilename );
    strcat( filename_pvd, ".pvd" );

    fpvtk = fopen( filename_pvd, "r" );

    while ( fgets( time_line, sizeof(time_line), fpvtk ) )
    {
      // Extract 8 first characters
      strncpy( start, time_line, 8 );

      // If 8 first characters equal "<DataSet", it is an output time line
      // We add to the vtk time series string
      if ( ! strcmp( start, start_ref ) )
        strcat( vtk_times_series, time_line );
    }

    fclose( fpvtk );
  }
}
*/


