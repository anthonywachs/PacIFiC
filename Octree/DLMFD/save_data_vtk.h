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
# ifndef PARAVIEW_VTU_MPIIO_WRITER 
#   define PARAVIEW_VTU_MPIIO_WRITER 0
# endif
# if PARAVIEW_VTU_MPIIO_WRITER && !_MPI
#   undef PARAVIEW_VTU_MPIIO_WRITER
#   define PARAVIEW_VTU_MPIIO_WRITER 0
# endif

# include "DLMFD_Output_vtu_foreach.h"
# include "vtkXMLHyperTreeGrid.h"


//----------------------------------------------------------------------------
void output_pvd( FILE* fp, char const* times_series )
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
void output_vtu_dlmfd_bndpts( RigidBody const* allrb, const int np,
	char const* fname )
//----------------------------------------------------------------------------
{
# if DLMFD_BOUNDARYPOINTS
    if ( pid() == 0 ) 
    {    
      int total_boundary_points = 0;
      for (size_t k = 0; k < np; k++) 
      {
        RigidBodyBoundary const* sbb = &(allrb[k].s);
	for (size_t j = 0; j < sbb->m; j++)
	  if ( sbb->deactivated[j] == 0 ) ++total_boundary_points;
      }     
  
      FILE* fdlm = fopen( fname, "w" ); 

      fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fdlm );
      fputs( "<UnstructuredGrid>\n", fdlm );
      fprintf( fdlm, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
      	total_boundary_points, total_boundary_points ); 
      fputs( "<Points>\n", fdlm );  
      fputs( "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
      	"format=\"ascii\">\n", fdlm );
      for (size_t k = 0; k < np; k++) 
      {
        RigidBodyBoundary const* sbb = &(allrb[k].s);
	int m = sbb->m;
	for (size_t j = 0; j < m; j++)
	  if ( sbb->deactivated[j] == 0 )
	  {  
	    fprintf( fdlm, "%g %g", sbb->x[j], sbb->y[j] );
#           if dimension == 3  
              fprintf( fdlm, " %g\n", sbb->z[j] );
#           else
              fprintf( fdlm, " 0.\n" );
#           endif	
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
      fputs( "<DataArray type=\"Int8\" Name=\"types\" "
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
void output_vtu_dlmfd_intpts( RigidBody const* allrb, const int np,
	char const* fname )
//----------------------------------------------------------------------------
{
# if DLMFD_INTERIORPOINTS
    int nb_ranks = 1;

#   if _MPI
      MPI_Comm_size( MPI_COMM_WORLD, &nb_ranks );
#   endif
     
    int number_interior_points = 0;
    int total_interior_points = 0;
    
    foreach(serial)
      if ( DLM_Flag[] < 1 && (int)DLM_Index.y[] > - 1 )
       	number_interior_points += 1;  
     	
    double* interior_coordx = (double*) calloc( number_interior_points, 
    	sizeof(double) ); 
    double* interior_coordy = (double*) calloc( number_interior_points, 
    	sizeof(double) );
    double* interior_coordz = (double*) calloc( number_interior_points, 
    	sizeof(double) );
    double* All_interior_coordx = NULL;
    double* All_interior_coordy = NULL;
    double* All_interior_coordz = NULL;
    int* interior_count = (int*) calloc( nb_ranks, sizeof(int) );    

    int counter = 0;
    foreach(serial)
      if ( DLM_Flag[] < 1 && (int)DLM_Index.y[] > - 1 )
      {
        interior_coordx[counter] = x;
        interior_coordy[counter] = y;
        interior_coordz[counter] = z;
        counter++;
      }
      
#   if _MPI
      int* interior_displace = (int*) calloc( nb_ranks, sizeof(int) );

      MPI_Allreduce( &number_interior_points, &total_interior_points, 1, 
      	MPI_INT, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allgather( &number_interior_points, 1, MPI_INT, interior_count, 1, 
    	MPI_INT, MPI_COMM_WORLD );
	
      int temp = 0;
      for(int i = 0; i < nb_ranks; i++)
      {
         interior_displace[i] = temp;
         temp += interior_count[i];
      }	
      
      if ( pid() == 0 ) 
      {
        All_interior_coordx = (double*) calloc( total_interior_points, 
		sizeof(double)); 
        All_interior_coordy = (double*) calloc( total_interior_points,
		sizeof(double));
        All_interior_coordz = (double*) calloc( total_interior_points, 
		sizeof(double));
      }
     
      MPI_Gatherv( interior_coordx, number_interior_points, MPI_DOUBLE, 
    	All_interior_coordx, interior_count, interior_displace, MPI_DOUBLE, 
	0,  MPI_COMM_WORLD );
      MPI_Gatherv( interior_coordy, number_interior_points, MPI_DOUBLE, 
    	All_interior_coordy, interior_count, interior_displace, MPI_DOUBLE, 
	0,  MPI_COMM_WORLD );
      MPI_Gatherv( interior_coordz, number_interior_points, MPI_DOUBLE, 
    	All_interior_coordz, interior_count, interior_displace, MPI_DOUBLE, 
	0,  MPI_COMM_WORLD ); 
	
      free( interior_displace ); 	     
#   else
      total_interior_points = number_interior_points;
      interior_count[0] = number_interior_points;
      All_interior_coordx = interior_coordx;
      All_interior_coordy = interior_coordy;
      All_interior_coordz = interior_coordz;       
#   endif
                     
    if ( pid() == 0 ) 
    {                    
      FILE* fdlm = fopen( fname, "w" ); 

      fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fdlm );
      fputs( "<UnstructuredGrid>\n", fdlm );
      fprintf( fdlm, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
      	total_interior_points, total_interior_points); 
      fputs( "<Points>\n", fdlm );  
      fputs( "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
      	"format=\"ascii\">\n", fdlm );
      	
      for(int i = 0; i < total_interior_points; i++)
      {
        fprintf( fdlm, "%g %g", All_interior_coordx[i], 
	  	All_interior_coordy[i] );
#       if dimension == 3  
          fprintf( fdlm, " %g\n", All_interior_coordz[i] );
#       else
          fprintf( fdlm, " 0.\n" );
#       endif	
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
      fputs( "<DataArray type=\"Int8\" Name=\"types\" "
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
    
    free( interior_coordx );
    free( interior_coordy );    
    free( interior_coordz );
    free( interior_count );        
#   if _MPI
      if ( pid() == 0 ) 
      {      
        free( All_interior_coordx );
        free( All_interior_coordy );      
        free( All_interior_coordz );
      }           
#   endif    
        
# endif     
}




//----------------------------------------------------------------------------
void save_data_vtk( scalar* list, vector* vlist, RigidBody const* allrb, 
	const int np, double const time )
//----------------------------------------------------------------------------
{
  static int cycle_number = 0; 
  if ( !cycle_number ) cycle_number = init_cycle_number;
  
  FILE * fpvtk; 
  char filename_pvd[80] = "";
  char suffix[80] = "";

# if PARAVIEW_VTU
    char filename_vtu[80] = "";            
    // Write the VTU file
    sprintf( filename_vtu, "%s", RESULT_DIR );
    strcat( filename_vtu, "/" );  
    strcat( filename_vtu, RESULT_FLUID_ROOTFILENAME );
#   if PARAVIEW_VTU_MPIIO_WRITER
      sprintf( suffix, "_T%d.vtu", cycle_number );
      strcat( filename_vtu, suffix );    
      if ( PARAVIEW_BINFILE ) 
        output_vtu_bin_foreach_MPIIO( list, vlist, filename_vtu );
      else output_vtu_ascii_foreach_MPIIO( list, vlist, filename_vtu );    
#   else   
#     if _MPI      
        sprintf( suffix, "_T%d_%d.vtu", cycle_number, pid() );
#     else
        sprintf( suffix, "_T%d.vtu", cycle_number );
#     endif	
      strcat( filename_vtu, suffix );
 
      fpvtk = fopen( filename_vtu, "w" );
      if ( PARAVIEW_BINFILE ) output_vtu_bin_foreach( list, vlist, fpvtk );
      else output_vtu_ascii_foreach( list, vlist, fpvtk );  
      fclose( fpvtk );
   
      // Write the PVTU file
#     if _MPI         
        if ( pid() == 0 ) 
        {
          char filename_pvtu[80] = ""; 
	  sprintf( filename_pvtu, "%s", RESULT_DIR );
          strcat( filename_pvtu, "/" );  
          strcat( filename_pvtu, RESULT_FLUID_ROOTFILENAME );    
          sprintf( suffix, "_T%d.pvtu", cycle_number );
          strcat( filename_pvtu, suffix );

          fpvtk = fopen( filename_pvtu, "w" );
    
          sprintf( filename_vtu, "%s", RESULT_FLUID_ROOTFILENAME );
          sprintf( suffix, "_T%d", cycle_number );
          strcat( filename_vtu, suffix );
          if ( PARAVIEW_BINFILE ) output_pvtu_bin( list, vlist, fpvtk, 
      		filename_vtu );
          else output_pvtu_ascii( list, vlist, fpvtk, filename_vtu );    

          fclose( fpvtk );
        }
#     endif	
#   endif
     
    // Write the PVD file  
    if ( pid() == 0 ) 
    {  
      sprintf( filename_pvd, "%s", RESULT_DIR );
      strcat( filename_pvd, "/" );  
      strcat( filename_pvd, RESULT_FLUID_ROOTFILENAME );
      strcat( filename_pvd, "_vtu.pvd" ); 

      fpvtk = fopen( filename_pvd, "w" );

      char time_line[200] = "";
      strcpy( time_line, "<DataSet timestep=" );
      sprintf( suffix, "\"%.4e\"", time );
      strcat( time_line, suffix );
      strcat( time_line, " group=\"\" part=\"0\" file=\"" );
      strcpy( filename_vtu, RESULT_FLUID_ROOTFILENAME );    
#     if PARAVIEW_VTU_MPIIO_WRITER || !_MPI
        sprintf( suffix, "_T%d.vtu", cycle_number );
#     else
        sprintf( suffix, "_T%d.pvtu", cycle_number );
#     endif      
      strcat( filename_vtu, suffix );
      strcat( time_line, filename_vtu );        
      strcat( time_line, "\"/>\n" );  
      strcat( vtu_field_times_series, time_line );    
      output_pvd( fpvtk, vtu_field_times_series );
  
      fclose( fpvtk );
    }
# endif


# if PARAVIEW_HTG 
    char filename_htg[80] = "";             
    // Write the HTG file
    sprintf( filename_htg, "%s", RESULT_DIR );
    strcat( filename_htg, "/" );  
    strcat( filename_htg, RESULT_FLUID_ROOTFILENAME );
    sprintf( suffix, "_T%d.htg", cycle_number );
    strcat( filename_htg, suffix );
    vtkXMLHyperTreeGrid *vtk_xml_hypertreegrid = NULL; 

#   if _MPI
      MPI_File fp;
      MPI_File_open( MPI_COMM_WORLD, filename_htg,
                MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp );
      vtk_xml_hypertreegrid =
      	vtk_xml_hypertreegrid_init( 2, PARAVIEW_DATATYPE_DOUBLE ? 9: 8, 1, 
		true, list, vlist, time );
      vtk_xml_hypertreegrid_to_file( vtk_xml_hypertreegrid, fp );
      vtk_xml_hypertreegrid_free( vtk_xml_hypertreegrid );
      MPI_File_close( &fp );
#   else
      FILE *fp = fopen( filename_htg, "w");
      if ( PARAVIEW_BINFILE )
        vtk_xml_hypertreegrid = vtk_xml_hypertreegrid_init( 2, 
		PARAVIEW_DATATYPE_DOUBLE ? 9: 8, 1, true, list, vlist, time );
      else 
        vtk_xml_hypertreegrid = vtk_xml_hypertreegrid_init( 2, 
		PARAVIEW_DATATYPE_DOUBLE ? 9: 8, 0, false, list, vlist, time );
      vtk_xml_hypertreegrid_to_file( vtk_xml_hypertreegrid, fp );
      fclose(fp);
      vtk_xml_hypertreegrid_free( vtk_xml_hypertreegrid );
#   endif
         
    // Write the PVD file  
    if ( pid() == 0 ) 
    {  
      sprintf( filename_pvd, "%s", RESULT_DIR );
      strcat( filename_pvd, "/" );  
      strcat( filename_pvd, RESULT_FLUID_ROOTFILENAME );
      strcat( filename_pvd, "_htg.pvd" ); 

      fpvtk = fopen( filename_pvd, "w" );

      char time_line[200] = "";
      strcpy( time_line, "<DataSet timestep=" );
      sprintf( suffix, "\"%.4e\"", time );
      strcat( time_line, suffix );
      strcat( time_line, " group=\"\" part=\"0\" file=\"" );
      strcpy( filename_htg, RESULT_FLUID_ROOTFILENAME );
      sprintf( suffix, "_T%d.htg", cycle_number );     
      strcat( filename_htg, suffix );        
      strcat( time_line, filename_htg );        
      strcat( time_line, "\"/>\n" );  
      strcat( htg_field_times_series, time_line );    
      output_pvd( fpvtk, htg_field_times_series );
  
      fclose( fpvtk );
    }
# endif
    
  
  // Write the last cycle number in a file for restart  
  if ( pid() == 0 ) 
  {
    char filename_lcn[256] = ""; 
    sprintf( filename_lcn, "%s", RESULT_DIR );
    strcat( filename_lcn, "/" );  
    strcat( filename_lcn, RESULT_FLUID_ROOTFILENAME );
    strcat( filename_lcn, "_lcn_vtk.txt" );

    fpvtk = fopen( filename_lcn, "w" );    

    fprintf( fpvtk, "%d\n", cycle_number );
    
    fclose( fpvtk );            
  }
  
# if PARAVIEW_DLMFD_BNDPTS
    char filename_bnd_vtu[80] = "";
    sprintf( filename_bnd_vtu, "%s", RESULT_DIR );
    strcat( filename_bnd_vtu, "/" );  
    strcat( filename_bnd_vtu, PARAVIEW_DLMFD_BNDPTS_FILENAME );  
    sprintf( suffix, "_T%d.vtu", cycle_number );
    strcat( filename_bnd_vtu, suffix );
        
    output_vtu_dlmfd_bndpts( allrb, np, filename_bnd_vtu );
  
    if ( pid() == 0 ) 
    {  
      sprintf( filename_pvd, "%s", RESULT_DIR );
      strcat( filename_pvd, "/" );  
      strcat( filename_pvd, PARAVIEW_DLMFD_BNDPTS_FILENAME );
      strcat( filename_pvd, ".pvd" ); 

      fpvtk = fopen( filename_pvd, "w" );

      char time_line[200] = "";
      strcpy( time_line, "<DataSet timestep=" );
      sprintf( suffix, "\"%.4e\"", time );
      strcat( time_line, suffix );
      strcat( time_line, " group=\"\" part=\"0\" file=\"" );
      strcat( time_line, filename_bnd_vtu );         
      strcat( time_line, "\"/>\n" );  
      strcat( vtk_bndpts_times_series, time_line );    
      output_pvd( fpvtk, vtk_bndpts_times_series );
  
      fclose( fpvtk );
    }         
# endif

# if PARAVIEW_DLMFD_INTPTS
    char filename_int_vtu[80] = "";
    sprintf( filename_int_vtu, "%s", RESULT_DIR );
    strcat( filename_int_vtu, "/" );  
    strcat( filename_int_vtu, PARAVIEW_DLMFD_INTPTS_FILENAME );  
    sprintf( suffix, "_T%d.vtu", cycle_number );
    strcat( filename_int_vtu, suffix );
        
    output_vtu_dlmfd_intpts( allrb, np, filename_int_vtu );
  
    if ( pid() == 0 ) 
    {  
      sprintf( filename_pvd, "%s", RESULT_DIR );
      strcat( filename_pvd, "/" );  
      strcat( filename_pvd, PARAVIEW_DLMFD_INTPTS_FILENAME );
      strcat( filename_pvd, ".pvd" ); 

      fpvtk = fopen( filename_pvd, "w" );

      char time_line[200] = "";
      strcpy( time_line, "<DataSet timestep=" );
      sprintf( suffix, "\"%.4e\"", time );
      strcat( time_line, suffix );
      strcat( time_line, " group=\"\" part=\"0\" file=\"" );
      strcat( time_line, filename_int_vtu );        
      strcat( time_line, "\"/>\n" );  
      strcat( vtk_intpts_times_series, time_line );    
      output_pvd( fpvtk, vtk_intpts_times_series );
  
      fclose( fpvtk );
    }         
# endif  
  
  ++cycle_number;       
}




//----------------------------------------------------------------------------
void reinitialize_vtk_restart( void )
//----------------------------------------------------------------------------
{
  // Get the last cycle cumber from previous simulation
  char filename_lcn[80] = "";
  sprintf( filename_lcn, "%s", RESULT_DIR );
  strcat( filename_lcn, "/" );  
  strcat( filename_lcn, RESULT_FLUID_ROOTFILENAME );
  strcat( filename_lcn, "_lcn_vtk.txt" );

  FILE * fpvtk = fopen( filename_lcn, "r" );    

  fscanf ( fpvtk, "%d", &init_cycle_number );
  ++init_cycle_number;
    
  fclose( fpvtk ); 
  
  // Re-initialize the time output series string
  if ( pid() == 0 ) 
  {    
    // Field files 
    char filename_pvd_root[80] = "";
    char time_line[256] = "";
    char start[9] = ""; 
    char start_ref[20] = "<DataSet";    
    sprintf( filename_pvd_root, "%s", RESULT_DIR );
    strcat( filename_pvd_root, "/" );  
    strcat( filename_pvd_root, RESULT_FLUID_ROOTFILENAME );
    
#   if PARAVIEW_VTU
      char filename_vtu_pvd[80] = "";
      strcpy( filename_vtu_pvd, filename_pvd_root );      
      strcat( filename_vtu_pvd, "_vtu.pvd" ); 

      fpvtk = fopen( filename_vtu_pvd, "r" ); 
    
      while ( fgets( time_line, sizeof(time_line), fpvtk ) ) 
      {      
        // Extract 8 first characters
        strncpy( start, time_line, 8 );
        start[8] = '\0';

        // If 8 first characters equal "<DataSet", it is an output time line
        // We add to the vtk time series string
        if ( ! strcmp( start, start_ref ) )
          strcat( vtu_field_times_series, time_line );      
      }
    
      fclose( fpvtk );
#   endif

#   if PARAVIEW_HTG
      char filename_htg_pvd[80] = "";
      strcpy( filename_htg_pvd, filename_pvd_root );      
      strcat( filename_htg_pvd, "_htg.pvd" ); 

      fpvtk = fopen( filename_htg_pvd, "r" ); 
    
      while ( fgets( time_line, sizeof(time_line), fpvtk ) ) 
      {      
        // Extract 8 first characters
        strncpy( start, time_line, 8 );
        start[8] = '\0';

        // If 8 first characters equal "<DataSet", it is an output time line
        // We add to the vtk time series string
        if ( ! strcmp( start, start_ref ) )
          strcat( htg_field_times_series, time_line );      
      }
    
      fclose( fpvtk );
#   endif          
    
#   if PARAVIEW_DLMFD_BNDPTS
      char filename_bnd_pvd[80] = "";
      sprintf( filename_bnd_pvd, "%s", RESULT_DIR );
      strcat( filename_bnd_pvd, "/" );  
      strcat( filename_bnd_pvd, PARAVIEW_DLMFD_BNDPTS_FILENAME );
      strcat( filename_bnd_pvd, ".pvd" ); 

      fpvtk = fopen( filename_bnd_pvd, "r" ); 
    
      while ( fgets( time_line, sizeof(time_line), fpvtk ) ) 
      {      
        // Extract 8 first characters
        strncpy( start, time_line, 8 );
        start[8] = '\0';

        // If 8 first characters equal "<DataSet", it is an output time line
        // We add to the vtk time series string
        if ( ! strcmp( start, start_ref ) )
          strcat( vtk_bndpts_times_series, time_line );      
      }
      
      fclose( fpvtk );             
#   endif 

#   if PARAVIEW_DLMFD_INTPTS
      char filename_int_pvd[80] = "";
      sprintf( filename_int_pvd, "%s", RESULT_DIR );
      strcat( filename_int_pvd, "/" );  
      strcat( filename_int_pvd, PARAVIEW_DLMFD_INTPTS_FILENAME );
      strcat( filename_int_pvd, ".pvd" ); 

      fpvtk = fopen( filename_int_pvd, "r" ); 
    
      while ( fgets( time_line, sizeof(time_line), fpvtk ) ) 
      {      
        // Extract 8 first characters
        strncpy( start, time_line, 8 );
        start[8] = '\0';

        // If 8 first characters equal "<DataSet", it is an output time line
        // We add to the vtk time series string
        if ( ! strcmp( start, start_ref ) )
          strcat( vtk_intpts_times_series, time_line );      
      }
      
      fclose( fpvtk );             
#   endif  
  }         
}
