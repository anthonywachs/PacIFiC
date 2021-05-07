/** 
# Wrapper for output functions with Paraview
*/
# include "dlmfd-output_vtu_foreach.h"


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
