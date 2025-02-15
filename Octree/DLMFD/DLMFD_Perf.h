/** 
# Performance of the Uzawa-DLMFD solver
*/

char dlmfd_perf_complete_name[80]; 

/** Outputs the performance of the DLMFD solver */
//----------------------------------------------------------------------------
void output_dlmfd_perf( timing* Uzawa, timing* Construction, const int i, 
	DLMFDptscells const* pp ) 
//----------------------------------------------------------------------------
{
  double mpitimings[npe()];
  
  static unsigned int iii = 0 ; 
  if ( iii == 0 )
  {  
    char buffer[80] = "";
    strcpy( buffer, RESULT_DIR );
    strcat( buffer, "/" );
    strcat( buffer, DLMFD_PERF_FILENAME );
    strcpy( dlmfd_perf_complete_name, buffer );
    ++iii;     
  }  
  
  static FILE* dlmfdperf;
  dlmfdperf = fopen( dlmfd_perf_complete_name, "a" );

  // Global timer/timings
  timing gns = timer_timing( perf.gt, i, perf.tnc, mpitimings );

  // Header
  fprintf( fout, "\n# DLMFD stats\n" );
  if ( pid() == 0 ) fprintf( dlmfdperf,"# DLMFD stats\n" );

  
  // DLMFD construction
  fprintf( fout, "#    Construction\n" );
  if ( pid() == 0 ) fprintf( dlmfdperf,"#    Construction\n" ); 
  
  // Number of cells
  fprintf( fout, "#       " GRIDNAME
	", %ld constrained cells, %ld constrained pts, %ld total cells\n", 
	pp->total_number_of_DLMFDcells/i, 
	pp->total_number_of_DLMFDpts/i, perf.tnc/i );
  if ( pid() == 0 ) fprintf( dlmfdperf, "#       " GRIDNAME
	", %ld constrained cells, %ld constrained pts, %ld total cells\n", 
	pp->total_number_of_DLMFDcells/i, 
	pp->total_number_of_DLMFDpts/i, perf.tnc/i );    
  
  // Display and write the statistics of the DLMFD Construction  
  fprintf( fout, "#       Time: %d steps, %g CPU, %.4g real\n",
	i, Construction->cpu, Construction->real );
  if ( pid() == 0 ) fprintf( dlmfdperf,"#       Time: %d steps, %g CPU, "
  	"%.4g real\n", i, Construction->cpu, Construction->real );
	
# if _MPI
    fprintf( fout,
	   "#       MPI: %d procs, min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   Construction->min, 100. * Construction->min / Construction->real,
	   Construction->avg, 100. * Construction->avg / Construction->real,
	   Construction->max, 100. * Construction->max / Construction->real );
    if ( pid() == 0 ) fprintf( dlmfdperf,
	   "#       MPI: %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   Construction->min, 100. * Construction->min / Construction->real,
	   Construction->avg, 100. * Construction->avg / Construction->real,
	   Construction->max, 100. * Construction->max / Construction->real );
# endif	   

  // Display and write the ratio DLMFD construction/total
  fprintf( fout, "#       Time to total time: " 
	"%g CPU, %.4g real\n", 
	Construction->cpu / gns.cpu, Construction->real / gns.real );
  if ( pid() == 0 ) fprintf( dlmfdperf, "#       Ratio to total: " 
	"%g CPU, %.4g real\n", 
	Construction->cpu / gns.cpu, Construction->real / gns.real );

# if _MPI  
    fprintf( fout,
	"#       MPI Time to total MPI time: min %.2g avg %.2g max %.2g\n",
	Construction->min / gns.min, Construction->avg / gns.avg, 
	Construction->max / gns.max );	
    if ( pid() == 0 ) fprintf( dlmfdperf,
	"#       MPI Time to total MPI time: min %.2g avg %.2g max %.2g\n",
	Construction->min / gns.min, Construction->avg / gns.avg, 
	Construction->max / gns.max );
# endif


  // DLMFD Uzawa
  fprintf( fout, "#    Uzawa\n" );
  if ( pid() == 0 ) fprintf( dlmfdperf,"#    Uzawa\n" ); 
    
  // Compute the speed with the total number of cells constrained by the DLMFD 
  // problem (that number is multiplied by the number of steps) 
  Uzawa->speed = pp->total_number_of_DLMFDcells / Uzawa->real;

  // Display and write the statistics of the DLMFD Uzawa solver  
  fprintf( fout,"#       Time: "
	"%d steps, %g CPU, %.4g real, %.3g points.step/s\n",
	i, Uzawa->cpu, Uzawa->real, Uzawa->speed );
  if ( pid() == 0 ) fprintf( dlmfdperf, "#       Time: "
	"%d steps, %g CPU, %.4g real, %.3g points.step/s\n",
	i, Uzawa->cpu, Uzawa->real, Uzawa->speed );

# if _MPI
    fprintf( fout,
	   "#       MPI: %d procs, min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   Uzawa->min, 100. * Uzawa->min / Uzawa->real,
	   Uzawa->avg, 100. * Uzawa->avg / Uzawa->real,
	   Uzawa->max, 100. * Uzawa->max / Uzawa->real );
    if ( pid() == 0 ) fprintf( dlmfdperf,
	   "#       MPI: %d procs, min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   Uzawa->min, 100. * Uzawa->min / Uzawa->real,
	   Uzawa->avg, 100. * Uzawa->avg / Uzawa->real,
 	   Uzawa->max, 100. * Uzawa->max / Uzawa->real );
# endif
  
  // Display and write the ratio DLMFD/total
  fprintf( fout, "#       Time to total time: "
	"%g CPU, %.4g real\n", 
	Uzawa->cpu / gns.cpu, Uzawa->real / gns.real );
  if ( pid() == 0 ) fprintf( dlmfdperf, "#       Time to total time: "
	"%g CPU, %.4g real\n", 
	Uzawa->cpu / gns.cpu, Uzawa->real / gns.real );
  
# if _MPI  
    fprintf( fout,
	"#       MPI Time to total MPI time: min %.2g "
	"avg %.2g max %.2g\n",
	Uzawa->min / gns.min, Uzawa->avg / gns.avg, Uzawa->max / gns.max );	
    if ( pid() == 0 ) fprintf( dlmfdperf,
	"#       MPI Time to total MPI time: min %.2g "
	"avg %.2g max %.2g\n",
	Uzawa->min / gns.min, Uzawa->avg / gns.avg, Uzawa->max / gns.max );
# endif


  // Construction and Uzawa together
  fprintf( fout, "#    All\n" );
  if ( pid() == 0 ) fprintf( dlmfdperf,"#    All\n" );  

  // Display and write the statistics of the DLMFD problem 
  fprintf( fout, "#       Time: "
	"%d steps, %g CPU, %.4g real\n",
	i, Uzawa->cpu + Construction->cpu, 
	Uzawa->real + Construction->real );
  if ( pid() == 0 ) fprintf( dlmfdperf, "#       Time: "
	"%d steps, %g CPU, %.4g real\n",
	i, Uzawa->cpu + Construction->cpu, 
	Uzawa->real + Construction->real );

# if _MPI
    fprintf( fout,
	   "#       MPI: %d procs, min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   Uzawa->min + Construction->min, 
	   100. * ( Uzawa->min + Construction->min ) / ( Uzawa->real 
	   	+ Construction->real ),
	   Uzawa->avg + Construction->avg, 
	   100. * ( Uzawa->avg + Construction->avg ) / ( Uzawa->real 
	   	+ Construction->real ),
	   Uzawa->max + Construction->max, 
	   100. * ( Uzawa->max + Construction->max ) / ( Uzawa->real 
	   	+ Construction->real ) );	
    if ( pid() == 0 ) fprintf( dlmfdperf,
	   "#       MPI: %d procs, min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   Uzawa->min + Construction->min, 
	   100. * ( Uzawa->min + Construction->min ) / ( Uzawa->real 
	   	+ Construction->real ),
	   Uzawa->avg + Construction->avg, 
	   100. * ( Uzawa->avg + Construction->avg ) / ( Uzawa->real 
	   	+ Construction->real ),
	   Uzawa->max + Construction->max, 
	   100. * ( Uzawa->max + Construction->max ) / ( Uzawa->real 
	   	+ Construction->real ) );
# endif
  
  // Display and write the ratio DLMFD/total
  fprintf( fout, "#       Time to total time: "
	"%g CPU, %.4g real\n", 
	( Uzawa->cpu + Construction->cpu ) / gns.cpu, 
	( Uzawa->real + Construction->real ) / gns.real );
  if ( pid() == 0 ) fprintf( dlmfdperf, "#       Time to total time: "
	"%g CPU, %.4g real\n", 
	( Uzawa->cpu + Construction->cpu ) / gns.cpu, 
	( Uzawa->real + Construction->real ) / gns.real );
  
# if _MPI  
    fprintf( fout,
	"#       MPI Time to total MPI time: min %.2g "
	"avg %.2g max %.2g\n",
	( Uzawa->min + Construction->min ) / gns.min, 
	( Uzawa->avg + Construction->avg ) / gns.avg, 
	( Uzawa->max + Construction->max ) / gns.max );	
    if ( pid() == 0 ) fprintf( dlmfdperf,
	"#       MPI Time to total MPI time: min %.2g "
	"avg %.2g max %.2g\n",
	( Uzawa->min + Construction->min ) / gns.min, 
	( Uzawa->avg + Construction->avg ) / gns.avg, 
	( Uzawa->max + Construction->max ) / gns.max );	
# endif
 
  // Write in the file with pointer dlmfdperf the total statistics
  if ( pid() == 0 ) fprintf( dlmfdperf,
	"\n# " GRIDNAME 
	", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
	i, gns.cpu, gns.real, gns.speed, (int) (datasize/sizeof(double)) );

# if _MPI
    if ( pid() == 0 ) fprintf( dlmfdperf,
	"# %d procs, MPI: min %.2g (%.2g%%) "
	"avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	npe(),
	gns.min, 100. * gns.min / gns.real,
	gns.avg, 100. * gns.avg / gns.real,
	gns.max, 100. * gns.max / gns.real );
# endif
  
  fflush( dlmfdperf ); 
}
