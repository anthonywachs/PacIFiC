/** 
# Wrapper for output functions with Paraview
*/
# include "DLMFD_Output_vtu_foreach.h"


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
void save_data( scalar * list, vector * vlist, double const time )
//----------------------------------------------------------------------------
{
  static int cycle_number = 0; 
  if ( !cycle_number ) cycle_number = init_cycle_number;
  
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




//----------------------------------------------------------------------------
void output_vtu_dlmfd_bpts( particle const* allpart, const int np,
	char const* fname )
//----------------------------------------------------------------------------
{
# if debugBD == 0
    if ( pid() == 0 ) 
    {
      char filename[80] = ""; 
      sprintf( filename, "%s", result_dir );
      strcat( filename, "/" );  
      strcat( filename, fname );
      strcat( filename, ".vtu" );
    
      int total_boundary_points = 0;
      for (int k = 0; k < np; k++) 
      {
        SolidBodyBoundary const* sbb = &(allpart[k].s);
        total_boundary_points += sbb->m;
      }     
  
      FILE* fdlm = fopen( filename, "w" ); 

      fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fdlm );
      fputs( "<UnstructuredGrid>\n", fdlm );
      fprintf( fdlm, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
      	total_boundary_points, total_boundary_points ); 
      fputs( "<Points>\n", fdlm );  
      fputs( "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
      	"format=\"ascii\">\n", fdlm );
      for (int k = 0; k < np; k++) 
      {
        SolidBodyBoundary const* sbb = &(allpart[k].s);
	int m = sbb->m;
	for (int j = 0; j < m; j++)
	{  
	  fprintf( fdlm, "%g %g", sbb->x[j], sbb->y[j] );
#         if dimension == 3  
            fprintf( fdlm, " %g\n", sbb->z[j] );
#         else
            fprintf( fdlm, " 0.\n" );
#         endif	
        }
      }
      fputs( "</DataArray>\n", fdlm );  
      fputs( "</Points>\n", fdlm );
      fputs( "<Cells>\n", fdlm );
      fputs( "<DataArray type=\"Int64\" Name=\"connectivity\" "
      	"format=\"ascii\">\n", fdlm );
      for (int j = 0; j < total_boundary_points; j++)
        fprintf( fdlm, "%d ", j );
      fprintf( fdlm, "\n" );
      fputs( "</DataArray>\n", fdlm ); 
      fputs( "<DataArray type=\"Int64\" Name=\"offsets\" "
      	"format=\"ascii\">\n", fdlm );
      for (int j = 0; j < total_boundary_points; j++)
        fprintf( fdlm, "%d ", j+1 );
      fprintf( fdlm, "\n" );		
      fputs( "</DataArray>\n", fdlm );
      fputs( "<DataArray type=\"Int64\" Name=\"types\" "
      	"format=\"ascii\">\n", fdlm );
      for (int j = 0; j < total_boundary_points; j++)
        fprintf( fdlm, "1 " );
      fprintf( fdlm, "\n" ); 	
      fputs( "</DataArray>\n", fdlm ); 
      fputs( "</Cells>\n", fdlm );
      fputs( "</Piece>\n", fdlm );
      fputs( "</UnstructuredGrid>\n", fdlm );            
      fputs( "</VTKFile>\n", fdlm );
                 	              	       
      fclose( fdlm );
    }
# endif     
}




 //----------------------------------------------------------------------------
void output_vtu_dlmfd_intpts( particle const* allpart, const int np,
	char const* fname )
//----------------------------------------------------------------------------
{
# if debugBD == 0
    int my_rank = 0, my_size = 1;

#   if _MPI
      MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
      MPI_Comm_size( MPI_COMM_WORLD, &my_size );
#   endif
     
    int number_interior_points = 0;
    int total_interior_points = 0;
    
    for (int k = 0; k < np; k++) 
      foreach()
        if ( flagfield[] < 1 && (int)index_lambda.y[] == k )
       	  number_interior_points += 1;
     
    double* interior_coordx = NULL;
    double* interior_coordy = NULL;
    double* interior_coordz = NULL;
    int* interior_count = NULL;
    int* interior_displace = NULL;
     	
    interior_coordx = (double*) calloc( number_interior_points, 
    	sizeof(double) ); 
    interior_coordy = (double*) calloc( number_interior_points, 
    	sizeof(double) );
    interior_coordz = (double*) calloc( number_interior_points, 
    	sizeof(double) );
    interior_count = (int*) calloc( my_size, sizeof(int) );
    interior_displace = (int*) calloc( my_size, sizeof(int) );

    int counter = 0;
    for (int k = 0; k < np; k++) 
      foreach()
     	if ( flagfield[] < 1 && (int)index_lambda.y[] == k ) 
       	{
          interior_coordx[counter] = x;
          interior_coordy[counter] = y;
          interior_coordz[counter] = z;
          counter++;
        }
      
    printf( "number_interior_points %d, counter is: %d, my rank is: %d\n", 
    	number_interior_points, counter, my_rank );
#   if _MPI
      MPI_Allreduce( &number_interior_points, &total_interior_points, 1, 
      	MPI_INT, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allgather( &number_interior_points, 1, MPI_INT, interior_count, 1, 
    	MPI_INT, MPI_COMM_WORLD );
#   else
      total_interior_points = number_interior_points;
      interior_count[0] = number_interior_points;
#   endif
    
    int temp = 0;
    for(int i = 0; i < my_size; i++)
    {
       interior_displace[i] = temp;
       temp += interior_count[i];
    }
    
    double* All_interior_coordx = NULL;
    double* All_interior_coordy = NULL;
    double* All_interior_coordz = NULL;
    
#   if _MPI
      All_interior_coordx = (double*) calloc( total_interior_points, 
    	sizeof(double)); 
      All_interior_coordy = (double*) calloc( total_interior_points,
    	sizeof(double));
      All_interior_coordz = (double*) calloc( total_interior_points, 
    	sizeof(double));
     
      MPI_Gatherv( interior_coordx, number_interior_points, MPI_DOUBLE, 
    	All_interior_coordx, interior_count, interior_displace, MPI_DOUBLE, 
	0,  MPI_COMM_WORLD );
      MPI_Gatherv( interior_coordy, number_interior_points, MPI_DOUBLE, 
    	All_interior_coordy, interior_count, interior_displace, MPI_DOUBLE, 
	0,  MPI_COMM_WORLD );
      MPI_Gatherv( interior_coordz, number_interior_points, MPI_DOUBLE, 
    	All_interior_coordz, interior_count, interior_displace, MPI_DOUBLE, 
	0,  MPI_COMM_WORLD );
#   else
      All_interior_coordx = interior_coordx;
      All_interior_coordy = interior_coordy;
      All_interior_coordz = interior_coordz;            
#   endif  
                     
    if ( pid() == 0 ) 
    {
      for(int i= 0; i < my_size; i++)
      {
   	printf( "Counts for each processor %d: %d\n", i, interior_count[i] );
   	printf( "Displacements for each processor %d: %d\n", i, 
		interior_displace[i] );
      }
    
      char filename[80] = ""; 
      sprintf( filename, "%s", result_dir );
      strcat( filename, "/" );  
      strcat( filename, fname );
      strcat( filename, ".vtu" );
                
      FILE* fdlm = fopen(filename, "w" ); 

      fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fdlm );
      fputs( "<UnstructuredGrid>\n", fdlm );
      fprintf( fdlm, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
      	total_interior_points, total_interior_points); 
      fputs( "<Points>\n", fdlm );  
      fputs( "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
      	"format=\"ascii\">\n", fdlm );
      	
      for (int k = 0; k < np; k++) 
      {
        printf( "Interior points are %d, and our counts are: %d", 
		allpart[k].Interior.n, total_interior_points );
	for(int i = 0; i < total_interior_points; i++)
     	{
          fprintf( fdlm, "%g %g", All_interior_coordx[i], 
	  	All_interior_coordy[i] );
#         if dimension == 3  
            fprintf( fdlm, " %g\n", All_interior_coordz[i] );
#         else
            fprintf( fdlm, " 0.\n" );
#         endif	
     	} 	  	
      }
      fputs( "</DataArray>\n", fdlm );  
      fputs( "</Points>\n", fdlm );
      fputs( "<Cells>\n", fdlm );
      fputs( "<DataArray type=\"Int64\" Name=\"connectivity\" "
      	"format=\"ascii\">\n", fdlm );
      for (int j = 0; j < total_interior_points; j++)
        fprintf( fdlm, "%d ", j );
      fprintf( fdlm, "\n" );
      fputs( "</DataArray>\n", fdlm ); 
      fputs( "<DataArray type=\"Int64\" Name=\"offsets\" "
      	"format=\"ascii\">\n", fdlm );
      for (int j = 0; j < total_interior_points; j++)
        fprintf( fdlm, "%d ", j+1 );
      fprintf( fdlm, "\n" );		
      fputs( "</DataArray>\n", fdlm );
      fputs( "<DataArray type=\"Int64\" Name=\"types\" "
      	"format=\"ascii\">\n", fdlm );
      for (int j = 0; j < total_interior_points; j++)
        fprintf( fdlm, "1 " );
      fprintf( fdlm, "\n" ); 	
      fputs( "</DataArray>\n", fdlm ); 
      fputs( "</Cells>\n", fdlm );
      fputs( "</Piece>\n", fdlm );
      fputs( "</UnstructuredGrid>\n", fdlm );            
      fputs( "</VTKFile>\n", fdlm );
                 	              	       
      fclose( fdlm );
    }
# endif     
}
