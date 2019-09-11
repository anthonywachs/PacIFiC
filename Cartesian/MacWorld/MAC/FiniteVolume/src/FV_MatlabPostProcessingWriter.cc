#include <MAC_Object.hh>
#include <FV_PostProcessingWriter.hh>
#include <FV_MatlabPostProcessingWriter.hh>
#include <FV_DomainAndFields.hh>
#include <FV_Mesh.hh>
#include <FV.hh>
#include <FV_DiscreteField.hh>
#include <MAC.hh>
#include <MAC_Bool.hh>
#include <MAC_Communicator.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>
#include <MAC_assertions.hh>
#include <FV_TimeIterator.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::cout ; 
using std::endl ;
using std::string ; 
using std::ostringstream ; 
using std::istringstream ; 
using std::ofstream ;
using std::ifstream ;
using std::ios ;


static size_t sizeof_Float32 = 4 ;
static size_t sizeof_Int32 = 4 ;

FV_MatlabPostProcessingWriter const* 
FV_MatlabPostProcessingWriter::PROTOTYPE = 
new FV_MatlabPostProcessingWriter( "matlab" ) ;


//----------------------------------------------------------------------
FV_MatlabPostProcessingWriter:: FV_MatlabPostProcessingWriter( 
      std::string const& a_name )
//----------------------------------------------------------------------
    : FV_PostProcessingWriter( a_name )
    , TIME_STRINGS( 0 )  
    , BUFFER( 0 )
    , ALLOCATED( 0 )
    , OFFSET( 0 )
    , WHOLE_EXTENT( 0 )
    , EXTENT( 0 )
    , EXTENT_AllProc( NULL )
    , BORDER( 0 )
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: FV_MatlabPostProcessingWriter" ) ;

}




//----------------------------------------------------------------------
FV_PostProcessingWriter*
FV_MatlabPostProcessingWriter:: create_replica( MAC_Object* a_owner, 
 		MAC_ModuleExplorer const* exp,
		MAC_Communicator const* com,
		list< FV_DiscreteField const* > a_fields,
		FV_Mesh const* a_primary_mesh,
		bool a_binary ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: create" ) ;

   FV_PostProcessingWriter* result 
   	= new FV_MatlabPostProcessingWriter( 
	  a_owner,  exp, com, a_fields, a_primary_mesh, a_binary ) ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_MatlabPostProcessingWriter:: FV_MatlabPostProcessingWriter(
		MAC_Object* a_owner, 
 		MAC_ModuleExplorer const* exp,
		MAC_Communicator const* com,
		list< FV_DiscreteField const* > a_fields,
		FV_Mesh const* a_primary_mesh,
		bool a_binary )
//----------------------------------------------------------------------
   : FV_PostProcessingWriter(a_owner)
   , EXP( exp )
   , COM( com )
   , FVFIELDS( a_fields )
   , PRIMARY_GRID( a_primary_mesh )
   , RES_DIRECTORY ( "Res" )
   , BASE_FILENAME( "save" )
   , TIME_FILENAME( "Res/save_time.xml" )
   , TIME_FILENAME_PVD( "Res/save.pvd" )
   , TIME_STRINGS( 0 )
   , CYCLE_NUMBER( 0 )
   , OUTPUT_DOMAIN_FILENAME("Res/current_domain.tmp")
   , OUTPUT_PROC_GRID_FILENAME("Res/current_proc_grid.tmp")
   , BINARY( a_binary )
   , BUFFER( 0 )
   , ALLOCATED( 0 )
   , OFFSET( 0 )
   , NB_SPACE_DIMENSION( a_primary_mesh->nb_space_dimensions() )
   , WHOLE_EXTENT( 0 )
   , EXTENT( 0 )
   , EXTENT_AllProc( NULL )
   , BORDER( 0 )
   , p_GLOBAL_MAX_INDEX( a_primary_mesh->get_global_max_index_in_domain() )
   , p_GLOBAL_MIN_INDEX( a_primary_mesh->get_global_min_index_in_domain() )
   , p_LOCAL_MAX_INDEX( 
   	a_primary_mesh->get_local_max_index_in_global_on_current_proc() )
   , p_LOCAL_MIN_INDEX( 
   	a_primary_mesh->get_local_min_index_in_global_on_current_proc() )
{
   MAC_LABEL( 
   "FV_MatlabPostProcessingWriter:: FV_MatlabPostProcessingWriter" ) ;
   
   //Global coordinates in domain
   BORDER.resize( 6 );
   BORDER(0) =  PRIMARY_GRID->get_main_domain_min_coordinate(0);
   BORDER(1) =  PRIMARY_GRID->get_main_domain_max_coordinate(0);
   BORDER(2) =  PRIMARY_GRID->get_main_domain_min_coordinate(1);
   BORDER(3) =  PRIMARY_GRID->get_main_domain_max_coordinate(1);
   if ( NB_SPACE_DIMENSION == 3 )
   {
     BORDER(4) = PRIMARY_GRID->get_main_domain_min_coordinate(2);
     BORDER(5) = PRIMARY_GRID->get_main_domain_max_coordinate(2);
   }
   

   // Global index in domain
   WHOLE_EXTENT.resize( 6 );
   WHOLE_EXTENT(0) = (*p_GLOBAL_MIN_INDEX)(0);
   WHOLE_EXTENT(1) = (*p_GLOBAL_MAX_INDEX)(0);
   WHOLE_EXTENT(2) = (*p_GLOBAL_MIN_INDEX)(1);   
   WHOLE_EXTENT(3) = (*p_GLOBAL_MAX_INDEX)(1);
   if ( NB_SPACE_DIMENSION == 3 )
   {
     WHOLE_EXTENT(4) = (*p_GLOBAL_MIN_INDEX)(2);   
     WHOLE_EXTENT(5) = (*p_GLOBAL_MAX_INDEX)(2);
   }
     
   // Local index in domain
   EXTENT.resize( 6 );
   EXTENT(0) = (*p_LOCAL_MIN_INDEX)(0);
   EXTENT(1) = (*p_LOCAL_MAX_INDEX)(0);
   EXTENT(2) = (*p_LOCAL_MIN_INDEX)(1);
   EXTENT(3) = (*p_LOCAL_MAX_INDEX)(1);
   if (NB_SPACE_DIMENSION == 3)
   {
     EXTENT(4) = (*p_LOCAL_MIN_INDEX)(2);
     EXTENT(5) = (*p_LOCAL_MAX_INDEX)(2);
   }

   NXPROCS = WHOLE_EXTENT(1) -  WHOLE_EXTENT(0) + 1;
   NYPROCS = WHOLE_EXTENT(3) -  WHOLE_EXTENT(2) + 1;
   NZPROCS = WHOLE_EXTENT(5) -  WHOLE_EXTENT(4) + 1;

   NXPROCS_LOC = PRIMARY_GRID->get_local_nb_points_on_current_proc( 0 );
   NYPROCS_LOC = PRIMARY_GRID->get_local_nb_points_on_current_proc( 1 );
   NZPROCS_LOC = 0 ;
   if ( NB_SPACE_DIMENSION == 3 )
     NZPROCS_LOC = PRIMARY_GRID->get_local_nb_points_on_current_proc( 2 );

   MY_COL =  PRIMARY_GRID->get_MPI_coordinates()[0];
   MY_ROW =  PRIMARY_GRID->get_MPI_coordinates()[1];
   MY_PLN =  0 ;
   if ( NB_SPACE_DIMENSION == 3 )
     MY_PLN = PRIMARY_GRID->get_MPI_coordinates()[2];

   // Get local min and max index of each processor
   size_t nb_ranks = COM->nb_ranks() ;
   EXTENT_AllProc = new vector< intVector >( nb_ranks, EXTENT ); 
   size_t iter = 0 ;
   intVector allExtent( nb_ranks * 6 ) ;
   COM->all_gather( EXTENT, allExtent ) ;

   for (size_t i=0;i<nb_ranks;++i)
     for( size_t j=0 ; j<6 ; j++, iter++ )
	 (*EXTENT_AllProc)[i](j) = allExtent(iter);

   // File names
   if ( exp->has_entry( "results_directory" ) )
     RES_DIRECTORY = exp->string_data( "results_directory" ); 
   if ( exp->has_entry( "files_rootname" ) )
     BASE_FILENAME = exp->string_data( "files_rootname" );         
   TIME_FILENAME = RES_DIRECTORY + "/" + BASE_FILENAME + "_time.xml" ;
   TIME_FILENAME_PVD = RES_DIRECTORY + "/" + "save.pvd" ;
   OUTPUT_DOMAIN_FILENAME = RES_DIRECTORY + "/" + "current_domain.tmp" ;
   OUTPUT_PROC_GRID_FILENAME =  RES_DIRECTORY + "/" + "current_proc_grid.tmp" ;

}




//----------------------------------------------------------------------
FV_MatlabPostProcessingWriter:: ~FV_MatlabPostProcessingWriter( void )
//----------------------------------------------------------------------
{
   MAC_LABEL(
   "FV_MatlabPostProcessingWriter:: ~FV_MatlabPostProcessingWriter" ) ;
   
   if( is_a_prototype() ) PROTOTYPE = 0 ;
//   EXTENT_AllProc->clear();
   delete EXTENT_AllProc;
   
}




//----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter::write_cycle( FV_TimeIterator const* t_it,
					size_t cycle_number )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: write_cycle" ) ;
   CYCLE_NUMBER = cycle_number;
   
   if ( COM->rank() == 0 )
   {
     // Save domain
     output_domain();

     // Save proc grid 
     output_proc_domain();
     
     // Time file
     std::ostringstream rootname ;
     rootname << BASE_FILENAME << "T" << CYCLE_NUMBER ;
     write_time_file( t_it, rootname.str() );
   }

   // Open files and test
   std::string fname = output_file_name( CYCLE_NUMBER, false, COM->rank() ) ;
   fname = RES_DIRECTORY + "/" + fname;
   std::ofstream file( fname.c_str(), std::ios::out ) ;
   if( !file )
   {
      MAC_Error::object()->raise_plain(
	 "unable to open the binary output file : " + fname ) ;
   }
   
   // 1st try : keep paraview format
   MAC_Module * matlab = MAC_Module::create( 0, "MATLAB" ) ;   
   build_module( matlab, false ) ;
   MAC_ModuleExplorer* mexp = MAC_ModuleExplorer::create( matlab, matlab ) ;   
   write_module( mexp, file, 0, false ) ; 
   file.close() ;
   matlab->destroy() ;
  
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: write_time_file(
	FV_TimeIterator const* t_it,
	std::string const& vtr_filename )
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: write_time_file" ) ;

   // First step : add a new line for new cycle if any
   std::string new_line = "<DataSet timestep=\"" ;
   std::ostringstream oss;
   oss << t_it->time() ; 

   new_line += oss.str() ;
   new_line += "\" file=\"" ;
   new_line += vtr_filename.c_str() ;
   new_line += "*\"/>" ;

   TIME_STRINGS.append( new_line ) ;
   
   // Second step : build time file with header, one line by cycle and end. 
   std::ofstream file( TIME_FILENAME.c_str(), std::ios::trunc ) ;
   if( !file )
   {
      MAC_Error::object()->raise_plain(
	 "unable to open the FTK output file : " + TIME_FILENAME ) ;
   }

   file << "<?xml version=\"1.0\"?>" << endl ;
   file << "</MATLABFile>" << endl << endl;
   file << "<FieldsOrdering>" << endl;
   list<FV_DiscreteField const*>::const_iterator il;
   for (il=FVFIELDS.begin(); il!=FVFIELDS.end() ; il++)
   {
     file << "<Field name=\"" << (*il)->name() << "\" nbcomps=\"" << 
     	(*il)->nb_components() << "\" location=\"" << (*il)->paraview_location()
	<< "\">" << endl;
   }  
   file << "</FieldsOrdering>" << endl << endl;   
   file << "<TimeSequence type=\"Collection\" version=\"0.1\""
        << " byte_order=\"LittleEndian\"" ;
   file << ">" << endl ;
   file << "<Collection>" << endl ;

   for( size_t i=0; i<TIME_STRINGS.size(); ++i )
   {
      file << TIME_STRINGS( i ) << endl ;
   }

   file << "</Collection>" << endl ;
   file << "</MATLABFile>" << endl ;
   file.close() ;
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: output_domain()
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: output_domain" ) ;
   
   std::ofstream file( OUTPUT_DOMAIN_FILENAME.c_str(), std::ios::trunc ) ;
   if( !file )
   {
      MAC_Error::object()->raise_plain(
	 "unable to open the output file : " + OUTPUT_DOMAIN_FILENAME ) ;
   }

   // Dimensions
   file << BORDER(0) << " " << BORDER(1) << " "
	<< BORDER(2) << " " << BORDER(3) << " "
	<< BORDER(4) << " " << BORDER(5) << " "  
	<< endl ;

   // Total number of points of the primary grid
   file << NXPROCS << " " << NYPROCS << " " <<  NZPROCS << " 0 0 0" << endl;
    
   // Primary grid coordinates by direction
   vector< doubleVector > const* maincoord = 
     	PRIMARY_GRID->get_global_main_coordinates();
   for (size_t i=0;i<NB_SPACE_DIMENSION;++i)
   {
     for (size_t j=WHOLE_EXTENT(2*i);j<=WHOLE_EXTENT(2*i+1);++j)
       file << (*maincoord)[i](j) << " " ;
     file << endl;  
   }

   // Total number of cell centres in the primary grid
   file << NXPROCS - 1 << " " << NYPROCS - 1 << " " << 
   	( NZPROCS == 1 ? 1 : NZPROCS - 1 ) << " 0 0 0" << endl;   
   
   // Primary grid cell centre coordinates by direction
   for (size_t i=0;i<NB_SPACE_DIMENSION;++i)
   {
     for (size_t j=WHOLE_EXTENT(2*i);j<WHOLE_EXTENT(2*i+1);++j)
       file << 0.5 * ( (*maincoord)[i](j) + (*maincoord)[i](j+1) ) << " " ;
     file << endl;  
   }       
             	
   file.close() ;
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: output_proc_domain()
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: output_domain" ) ;
   
   std::ofstream file( OUTPUT_PROC_GRID_FILENAME.c_str(), std::ios::trunc ) ;
   if( !file )
   {
      MAC_Error::object()->raise_plain(
	 "unable to open the output file : " + OUTPUT_PROC_GRID_FILENAME) ;
   }

   file << PRIMARY_GRID->get_domain_decomposition()[0] << " "
   	<< PRIMARY_GRID->get_domain_decomposition()[1] << " "
   	<< ( NB_SPACE_DIMENSION == 3 ?  
		PRIMARY_GRID->get_domain_decomposition()[2] : 0 ) << " "
   	<< 0 << " " << 0 << " " << 0  
   	<< endl;

   for (size_t i=0;i<COM->nb_ranks();++i)
   {
     for( size_t j=0 ; j<6 ; j++) file << (*EXTENT_AllProc)[i](j) << " " ; 
     file << endl ;
   }
   
   // Primary grid coordinates by direction on each proc
   vector< doubleVector > const* maincoord = 
     	PRIMARY_GRID->get_global_main_coordinates();
   for (size_t i=0;i<COM->nb_ranks();++i)   
     for (size_t m=0;m<NB_SPACE_DIMENSION;++m)
     {
       for (int j=(*EXTENT_AllProc)[i](2*m);j<=(*EXTENT_AllProc)[i](2*m+1);++j)
         file << (*maincoord)[m](j) << " " ;
       file << endl;  
     }           
   
   for (size_t i=0;i<COM->nb_ranks();++i)
   {
     for( size_t j=0 ; j<3 ; j++) 
       file << (*EXTENT_AllProc)[i](2*j) << " " << 
       	( (*EXTENT_AllProc)[i](2*j+1) == 0 ? 0 : 
		(*EXTENT_AllProc)[i](2*j+1) - 1 ) << " "; 
     file << endl ;
   }
   
   // Primary grid cell centre coordinates by direction on each proc
   for (size_t i=0;i<COM->nb_ranks();++i)   
     for (size_t m=0;m<NB_SPACE_DIMENSION;++m)
     {
       for (int j=(*EXTENT_AllProc)[i](2*m);j<(*EXTENT_AllProc)[i](2*m+1);++j)
         file << 0.5 * ( (*maincoord)[m](j) + (*maincoord)[m](j+1) ) << " " ;
       file << endl;  
     }   

   file.close() ;
}




//-----------------------------------------------------------------------
std::string
FV_MatlabPostProcessingWriter:: output_file_name( size_t nb,
                                  bool parallel,
                                  size_t rank )
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: output_file_name" ) ;
   MAC_CHECK( nb<99999 ) ;
   std::string file_extension ;
   if ( BINARY )
     file_extension = ( parallel ? ".pbin" : ".bin" ) ;
   else
     file_extension = ( parallel ? ".pm" : ".m" ) ; 
   
   std::ostringstream tmp ;
   tmp << BASE_FILENAME ;
   tmp << "T" << nb ;
   tmp << "_" << rank ;
   
   std::string result = tmp.str() + file_extension ;

   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: build_module(
				MAC_Module* matlab,
                         	bool parallel ) 
//----------------------------------------------------------------------
{
   MAC_LABEL("FV_MatlabPostProcessingWriter:: build_module") ;

   std::string pre = ( parallel ? "P" : "" ) ;
   MAC_Module* grid = MAC_Module::create( matlab, pre+"RectilinearGrid" ) ;
   MAC_Module* piece = MAC_Module::create( grid, "Piece" ) ;
   MAC_Module* base = ( parallel ? grid : piece ) ;
      
   // Create PointData and CellData modules
   MAC_Module* PointData = MAC_Module::create( base, pre+"PointData" ) ;
   MAC_Module* CellData = MAC_Module::create( base, pre+"CellData" ) ;
   
   // Write fields in files
   list<FV_DiscreteField const*>::const_iterator il;
   for (il=FVFIELDS.begin(); il!=FVFIELDS.end() ; il++)
       (*il)->write_field( PointData, CellData );
   
   base->add_module( PointData ) ;
   base->add_module( CellData ) ;   
   
   if( !parallel ) grid->add_module(piece) ;
   matlab->add_module(grid) ;   
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: write_module( MAC_ModuleExplorer* matlab,
                           std::ofstream& file,
                           size_t level,
                           bool parallel )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: write_module" ) ;

     
   for( matlab->start_entry_iterator() ;
        matlab->is_valid_entry() ;
        matlab->go_next_entry() ) 
     {
       MAC_Data* data = matlab->data( 0 ) ;
       MAC_Data::Type dt = data->data_type() ;
       if( dt==MAC_Data::DoubleArray2D || dt==MAC_Data::DoubleVector ) 
	 {
	if( !parallel )
         {           
            if( dt==MAC_Data::DoubleArray2D )
            {
               doubleArray2D const& tab = data->to_double_array2D() ;
               size_t nbc = data->to_double_array2D().index_bound(0) ;
               size_t ncmax = ( nbc==1 ? 1 : 3 ) ;
	       
	       //Warning: change Paraview-index 'ij' to Matlab index 'ji' 
	       for( size_t j=0 ; j<ncmax ; j++ )
	       {
		 //Here format="Unformatted" Matlab cf. notice
		 //+ add some header parameters s.a. nxprocs,nyprocs,nzprocs 
		 start_output( sizeof_Float32, ncmax*tab.index_bound(1) ) ;
		 for( size_t i=0 ; i<tab.index_bound(1) ; i++ )
		   write_double( file, ( j<nbc ? tab(j,i) : 0.0 ) ) ;
	       }
	    }
            else if( dt==MAC_Data::DoubleVector  )
            {
               doubleVector const& tab = data->to_double_vector() ;
	       if (tab.size() > 0)
	       {
		 start_output( sizeof_Float32, tab.size() ) ;
	         for( size_t i=0 ; i<tab.size() ; i++ )
	           write_double( file, tab(i) ) ;
	       }
            }
         }
         
      }
      else if( ! ( dt==MAC_Data::String || dt==MAC_Data::Int
      				|| dt==MAC_Data::IntVector) ) 
      {
         MAC_Error::object()->raise_internal(
            " Bad type "+data->data_type() ) ;
      }
      data->destroy() ;
   }

   for( matlab->start_module_iterator() ;
        matlab->is_valid_module() ;
        matlab->go_next_module() ) 
   {
      MAC_ModuleExplorer* sexp = matlab->create_subexplorer( 0 ) ;
      write_module( sexp, file, level+1, parallel ) ;
      sexp->destroy() ;
   }
   
   if( level == 0 && BINARY && !parallel ) 
   {
     file.write( BUFFER, OFFSET ) ;
     delete [] BUFFER ; BUFFER = 0 ;
     ALLOCATED = 0 ;
     OFFSET = 0 ;
   }
   
} 




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: start_output( size_t size,
							size_t number )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: start_output" ) ;
   if( BINARY )
   {
     int current_output_size = size*number ;
     unsigned long ncomp = current_output_size 
       + (current_output_size+999)/1000 + 12 + sizeof_Int32 ;
     check_allocated( ncomp ) ;
     header_output();
     
   }
   
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: header_output( void )  
//-----------------------------------------------------------------------
{
  MAC_LABEL( "FV_MatlabPostProcessingWriter:: header_output" ) ;

  if( BINARY )
  {
    float rank;
    rank = float(COM->rank());
    
    // Rank id
    store_double(rank);
    
    //Resolution by proc domain
    store_double(float(NXPROCS_LOC));
    store_double(float(NYPROCS_LOC));
    if ( NB_SPACE_DIMENSION == 3 )
      store_double(float(NZPROCS_LOC));
    
    //MPI Cart.
    store_double(float(MY_COL));
    store_double(float(MY_ROW));
    if ( NB_SPACE_DIMENSION == 3 )
      store_double(float(MY_PLN));
  }
  
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: write_double( std::ofstream& file,
							double val )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: write_double" ) ;
   MAC_CHECK( sizeof_Float32==sizeof(float) ) ;
   
   if( BINARY )
   {
      MAC_CHECK( OFFSET + sizeof_Float32 <= ALLOCATED ) ;
      *((float*)&(BUFFER[OFFSET])) = (float)val ;
      OFFSET += sizeof_Float32  ;
   }
   else
      file << " " << (float)val ;
}




//-----------------------------------------------------------------------
size_t
FV_MatlabPostProcessingWriter:: store_double( double val )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: store_double" ) ;
   MAC_CHECK( OFFSET + sizeof_Float32 <= ALLOCATED ) ;
   
   size_t result = OFFSET ;
   *((float*)&(BUFFER[OFFSET])) = val ;
   OFFSET += sizeof_Float32  ;
   return result ;
}





//-----------------------------------------------------------------------
size_t
FV_MatlabPostProcessingWriter:: store_int( int val )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: store_int" ) ;
   MAC_CHECK( OFFSET + sizeof_Int32 <= ALLOCATED ) ;
   
   size_t result = OFFSET ;
   *((int*)&(BUFFER[OFFSET])) = val ;
   OFFSET += sizeof_Int32  ;
   return result ;
}



//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: check_allocated( size_t size )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: check_allocated" ) ;
   if( OFFSET + size >= ALLOCATED ) 
   {
      size_t new_size = MAC::max( 2*ALLOCATED, (size_t)1024 ) ;
      new_size = MAC::max( new_size, 2*(OFFSET+size) ) ;
      new_size = 4 * ( new_size/4 +1 ) ; // allignement sur 4 bytes
      
      char * new_buffer = new char [ new_size ] ;
      for( size_t i=0 ;i<OFFSET ;i++ )
         new_buffer[i] = BUFFER[i] ;
      if( BUFFER!=0 ) delete [] BUFFER ;
      BUFFER = new_buffer ;
      ALLOCATED = new_size ;
      
   }
   MAC_CHECK_POST( OFFSET+size<ALLOCATED ) ;
}




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: clearResultFiles( void )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: clearResultFiles" ) ;

   if ( MAC_Exec::communicator()->rank() == 0 )
   {
     FV::out() << "   CLEAR MATLAB RESULTS FILES IN DIRECTORY " << 
     	RES_DIRECTORY << endl;

     char* path = getenv("MAC_HOME");
     string exe(path);
     exe = exe+"/tools/ShellScripts/FVMatlab_clear_exec "+RES_DIRECTORY;
     system(exe.c_str());
   }
}




//-----------------------------------------------------------------------
size_t
FV_MatlabPostProcessingWriter:: getPreviousCycleNumber( void )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: getPreviousCycleNumber" ) ; 
   
   string tline,previous_tline,buffer,part;

   if ( MAC_Exec::communicator()->rank() == 0 )  
      FV::out() << "   GET MATLAB PREVIOUS CYCLE NUMBER ";
   
   //if ( PostProcessor == "matlab" )
   //{
     // Read pvd file
   ifstream fileIN( TIME_FILENAME.c_str(), ios::in );
   if(fileIN) //if fileIN exist
   {
     while (tline != "</Collection>")
     {
       previous_tline = tline;
       getline( fileIN, tline );
     }
     // Extract cycle number from last result line in pvd file
     istringstream iss( previous_tline ); 
     iss >> buffer >> buffer >> part;
     size_t pos = part.find( "T" );
     string sub = part.substr( pos );
     sub.erase( sub.begin(), sub.begin()+1 );
     // Remark: total extension in .xml file is not specified after the cycle number
     // and hence is written as "*"/>", i.e. 4 characters
     size_t ncerase = 4 ;   
     sub.erase( sub.end()-ncerase, sub.end() ); 
     istringstream issNum( sub );
     issNum >> CYCLE_NUMBER;
   }
   //}
   //else

   //{
   //}
   else
   {
     // Read pvd file
     ifstream fileIN( TIME_FILENAME_PVD.c_str(), ios::in );
     while (tline != "</Collection>")
     {
       previous_tline = tline;
       getline( fileIN, tline );      
     }
     // Extract cycle number from last result line in pvd file
     istringstream iss( previous_tline ); 
     iss >> buffer >> buffer >> buffer >> buffer >> part;
     size_t pos = part.find( "T" );
     string sub = part.substr( pos );
     sub.erase( sub.begin(), sub.begin()+1 );
     // Remark: 
     // if serial, extension of results file is ".vtr", hence total extension
     // in .pvd file is ".vtr"/>", i.e. 7 characters
     // if parallel, extension of results file is ".pvtr", hence total extension
     // in .pvd file is ".pvtr"/>", i.e. 8 characters
     size_t ncerase = MAC_Exec::communicator()->nb_ranks() == 1 ? 7 : 8 ;   
     sub.erase( sub.end()-ncerase, sub.end() ); 
     istringstream issNum( sub );
     issNum >> CYCLE_NUMBER;

   }
   //}

   if ( MAC_Exec::communicator()->rank() == 0 )  
     FV::out() << CYCLE_NUMBER << endl;  
   return CYCLE_NUMBER ;    
}  




//-----------------------------------------------------------------------
void
FV_MatlabPostProcessingWriter:: readTimeFile( FV_TimeIterator const* t_it,
      			size_t& cycle_number  )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_MatlabPostProcessingWriter:: readTimeFile" ) ; 

   string tline;
   bool b_startcollection = false ;
   stringVector TIME_STRINGS_temp( 1 ) ;
   size_t count = 0 ;
   
   ifstream fileIN( TIME_FILENAME.c_str(), ios::in );
   // Switch for restart simulation with Matlab postprocessing
   if(fileIN) // if fileIN exist
   {
     getline( fileIN, tline );
     while ( !fileIN.eof() ) 
     { 
       if ( tline != "</Collection>" && tline != "</MATLABFile>" 
     	  && b_startcollection ) 
       {
	 if ( count ) TIME_STRINGS_temp.append( tline ) ;
	 else TIME_STRINGS_temp(0) = tline ;
	 ++count ;
       }
       if ( tline == "<Collection>" ) b_startcollection = true ;
       getline( fileIN, tline );    
     }
     fileIN.close(); 
   
     // Append all results lines but last, since it is written at restart
     // This avoids to have the last line twice in the pvd file
     for (size_t i=0; i<TIME_STRINGS_temp.size()-1;++i)
       TIME_STRINGS.append( TIME_STRINGS_temp(i) );    
   }
   else
   {
     string tline, keyword ;
     double starting_time = t_it->time(), time_step = t_it->time_step(),
   	  output_time ; 
     list<double> output_times ;
     list<string> output_time_lines ;
     ifstream fileIN( TIME_FILENAME_PVD.c_str(), ios::in );
     getline( fileIN, tline );
     while ( !fileIN.eof() ) 
     { 
       istringstream iss( tline );
       iss >> keyword;
       if ( keyword == "<DataSet" )
       {
         iss >> keyword;
         // delete the last "
         keyword.erase( keyword.end()-1, keyword.end() );
         // delete the timestep=" 
         keyword.erase( keyword.begin(), keyword.begin()+10 );
         // Convert to double
         istringstream issd( keyword );
         issd >> output_time ;
         output_times.push_back( output_time );
         output_time_lines.push_back( tline );	  
       }  
       getline( fileIN, tline );    
     }
   fileIN.close(); 
   
   // Test ultimate and penultimate output times vs starting_time
   list<double>::iterator id = output_times.end() ;
   list<string>::iterator il = output_time_lines.end() ;   
   --id; --il;
   if ( fabs( *id - starting_time ) < time_step )
   {
     if ( COM->rank() == 0 )
       FV::out() << "   Starting time matches last output time in save.pvd"
       		<< endl;
   }
   else if ( output_times.size() > 1 )
   {
      --id; --il;
      if ( fabs( *id - starting_time ) < time_step )
      {
        if ( COM->rank() == 0 )
	{
          FV::out() << "   Starting time matches penultimate output time in "
	  	"save.pvd" << endl;
	  --cycle_number;
          FV::out() << "   PARAVIEW PREVIOUS CYCLE NUMBER IS " << cycle_number
	  	<< endl;
	}
	output_time_lines.pop_back();
      }
      else if ( COM->rank() == 0 )
        FV::out() << "   WARNING : Starting time does match any of the "
		"previous output times\n   (neither ultimate nor penultimate)" 
		<< endl;        
   }
   else if ( COM->rank() == 0 )
     FV::out() << "   WARNING : Starting time does match any of the "
		"previous output times\n   (neither ultimate nor penultimate)" 
		<< endl;

   
   // Append all results lines but last, since it is written at restart
   // This avoids to have the last line twice in the pvd file
   output_time_lines.pop_back();
   for (il=output_time_lines.begin();il!=output_time_lines.end();il++)
     TIME_STRINGS.append( *il );      
   }
}   
