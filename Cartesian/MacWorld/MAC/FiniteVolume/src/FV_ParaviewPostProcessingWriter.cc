#include <MAC_Object.hh>
#include <FV_ParaviewPostProcessingWriter.hh>
#include <FV_PostProcessingWriter.hh>
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
#include <zlib.h>
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


FV_ParaviewPostProcessingWriter const* 
FV_ParaviewPostProcessingWriter::PROTOTYPE = 
new FV_ParaviewPostProcessingWriter( "paraview" ) ;


//----------------------------------------------------------------------
FV_ParaviewPostProcessingWriter:: FV_ParaviewPostProcessingWriter( 
         std::string const& a_name )
//----------------------------------------------------------------------
    : FV_PostProcessingWriter( a_name )
    , PVD_STRINGS( 0 )  
    , BUFFER( 0 )
    , ALLOCATED( 0 )
    , OFFSET( 0 )
    , WHOLE_EXTENT( 0 )
    , EXTENT( 0 )
    , EXTENT_AllProc( NULL )
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: " 
   	"FV_ParaviewPostProcessingWriter" ) ;

}

//----------------------------------------------------------------------
FV_PostProcessingWriter*
FV_ParaviewPostProcessingWriter:: create_replica( MAC_Object* a_owner, 
 		MAC_ModuleExplorer const* exp,
		MAC_Communicator const* com,
		list< FV_DiscreteField const* > a_fields,
		FV_Mesh const* a_primary_mesh,
		bool a_binary ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: create_replica" ) ;
    
   FV_PostProcessingWriter* result
   	= new FV_ParaviewPostProcessingWriter( 
	       a_owner, exp, com, a_fields, a_primary_mesh, a_binary ) ;

   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_ParaviewPostProcessingWriter:: FV_ParaviewPostProcessingWriter(
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
   , PVD_FILENAME( "Res/save.pvd" )
   , PVD_STRINGS( 0 )
   , CYCLE_NUMBER( 0 )
   , BINARY( a_binary )
   , BUFFER( 0 )
   , ALLOCATED( 0 )
   , OFFSET( 0 )
   , NB_SPACE_DIMENSION( a_primary_mesh->nb_space_dimensions() )
   , WHOLE_EXTENT( 0 )
   , EXTENT( 0 )
   , EXTENT_AllProc( NULL )
   , p_GLOBAL_MAX_INDEX( a_primary_mesh->get_global_max_index_in_domain() )
   , p_GLOBAL_MIN_INDEX( a_primary_mesh->get_global_min_index_in_domain() )
   , p_LOCAL_MAX_INDEX( 
   	a_primary_mesh->get_local_max_index_in_global_on_current_proc() )
   , p_LOCAL_MIN_INDEX( 
   	a_primary_mesh->get_local_min_index_in_global_on_current_proc() )
{
   MAC_LABEL( 
   "FV_ParaviewPostProcessingWriter:: FV_ParaviewPostProcessingWriter" ) ;
   
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
     
   EXTENT.resize( 6 );
   EXTENT(0) = (*p_LOCAL_MIN_INDEX)(0);
   EXTENT(1) = (*p_LOCAL_MAX_INDEX)(0);
   EXTENT(2) = (*p_LOCAL_MIN_INDEX)(1);
   EXTENT(3) = (*p_LOCAL_MAX_INDEX)(1);
   if ( NB_SPACE_DIMENSION == 3 )
   {
     EXTENT(4) = (*p_LOCAL_MIN_INDEX)(2);
     EXTENT(5) = (*p_LOCAL_MAX_INDEX)(2);
   }

   // Get local min and max index of each processor
   EXTENT_AllProc = new vector< intVector >(COM->nb_ranks(),EXTENT); 
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t nb_ranks = macCOMM->nb_ranks(), iter=0;
   intVector allExtent( nb_ranks*6 ) ;
   macCOMM->all_gather( EXTENT, allExtent ) ;

   for (size_t i=0;i<COM->nb_ranks();++i)
     for( size_t j=0 ; j<6 ; j++, iter++ )
       (*EXTENT_AllProc)[i](j) = allExtent(iter);
       
   // File names
   if ( exp->has_entry( "results_directory" ) )
     RES_DIRECTORY = exp->string_data( "results_directory" ); 
   if ( exp->has_entry( "files_rootname" ) )
     BASE_FILENAME = exp->string_data( "files_rootname" );         
   PVD_FILENAME = RES_DIRECTORY + "/" + BASE_FILENAME + ".pvd" ;
}




//----------------------------------------------------------------------
FV_ParaviewPostProcessingWriter:: ~FV_ParaviewPostProcessingWriter( void )
//----------------------------------------------------------------------
{
   MAC_LABEL(
   "FV_ParaviewPostProcessingWriter:: ~FV_ParaviewPostProcessingWriter" ) ;
   
 if( is_a_prototype() ) PROTOTYPE = 0 ;
 //EXTENT_AllProc->clear();
 delete EXTENT_AllProc;

}




//----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter::write_cycle( FV_TimeIterator const* t_it,
					size_t cycle_number )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: write_cycle" ) ;
   CYCLE_NUMBER = cycle_number;

   // Write PVD file
   if( COM->nb_ranks()>1 && COM->rank()==0 )
   {
      write_pvd_file( t_it, output_file_name( CYCLE_NUMBER, true, 0 ) ) ;
   }
   else if( COM->rank()==0 )
   {
      write_pvd_file( t_it, output_file_name( CYCLE_NUMBER, false, 0 ) ) ;
   }

   // Open vtr files and test
   std::string fname = output_file_name( CYCLE_NUMBER, false, COM->rank() ) ;
   fname = RES_DIRECTORY + "/" + fname;
   std::ofstream file( fname.c_str(), std::ios::out ) ;
   if( !file )
   {
      MAC_Error::object()->raise_plain(
	 "unable to open the VTR output file : " + fname ) ;
   }
   
   MAC_Module * vtk = MAC_Module::create( 0, "VTKFile" ) ;   
   build_vtr( vtk, false ) ;
   MAC_ModuleExplorer* mexp = MAC_ModuleExplorer::create( vtk, vtk ) ;   
   write_vtk( mexp, file, 0, false ) ;                                
   file.close() ;
   vtk->destroy() ;

   if( COM->nb_ranks() > 1 && COM->rank()==0 ) 
   {
      vtk = MAC_Module::create( 0, "VTKFile" ) ;
      build_vtr( vtk, true ) ;
      MAC_Module* grid = vtk->module( "PRectilinearGrid" ) ;
      grid->add_entry( "GhostLevel", MAC_Int::create( grid, 0 ) ) ;
      for( size_t i=0 ; i<COM->nb_ranks() ; i++ ) 
      {
         std::ostringstream npiece ;
         npiece << "Piece#" << i ;
         MAC_Module* piece = MAC_Module::create( grid, npiece.str() ) ;
         std::string file_name = output_file_name( CYCLE_NUMBER, false, i ) ;
         piece->add_entry( "Extent",
	 	MAC_IntVector::create( piece, (*EXTENT_AllProc)[i] ) ) ;
         piece->add_entry( "Source", MAC_String::create( piece, file_name ) ) ;
         grid->add_module( piece ) ;
      }
      
      // Open and write .pvtr files"
      std::string pfname = RES_DIRECTORY + "/" +
		output_file_name( CYCLE_NUMBER, true, 0 ) ;
      
      std::ofstream pfile( pfname.c_str() ) ;
      if( !pfile )
      {
         MAC_Error::object()->raise_plain(
            "unable to open the VTK output file : " + pfname ) ;
      }
      mexp = MAC_ModuleExplorer::create( vtk, vtk ) ;
      write_vtk( mexp, pfile, 0, true ) ;
      pfile.close() ;
      vtk->destroy() ;
   }
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: write_pvd_file(
				FV_TimeIterator const* t_it,
                                std::string const& vtr_filename )
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: write_pvd_file" ) ;

   // First step : add a new line for new cycle if any
   std::string new_line = "<DataSet timestep=\"" ;
   std::ostringstream oss;
   oss << t_it->time() ; 

   new_line += oss.str() ;
   new_line += "\" group=\"\" part=\"0\" file=\"" ;
   new_line += vtr_filename.c_str() ;
   new_line += "\"/>" ;

   PVD_STRINGS.append( new_line ) ;
   
   // Second step : build pvd file with header, one line by cycle and end. 
   std::ofstream file( PVD_FILENAME.c_str(), std::ios::trunc ) ;
   if( !file )
   {
      MAC_Error::object()->raise_plain(
	 "unable to open the FTK output file : " + PVD_FILENAME ) ;
   }

   file << "<?xml version=\"1.0\"?>" << endl ;
   file << "<VTKFile type=\"Collection\" version=\"0.1\""
        << " byte_order=\"LittleEndian\"" ;
   file << ">" << endl ;
   file << "<Collection>" << endl ;

   for( size_t i=0; i<PVD_STRINGS.size(); ++i )
   {
      file << PVD_STRINGS( i ) << endl ;
   }

   file << "</Collection>" << endl ;
   file << "</VTKFile>" << endl ;
   file.close() ;
}




//-----------------------------------------------------------------------
std::string
FV_ParaviewPostProcessingWriter:: output_file_name( size_t nb,
                                  bool parallel,
                                  size_t rank )
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: output_file_name" ) ;
   MAC_CHECK( nb<99999 ) ;
   std::string file_extension = ( parallel ? ".pvtr" : ".vtr" ) ;
   
   std::ostringstream tmp ;
   tmp << BASE_FILENAME ;
   tmp << "T" << nb ;
   if( !parallel && COM->nb_ranks() > 1 )
   {
      tmp << "_" << rank ;
   }
   
   std::string result = tmp.str() + file_extension ;

   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: build_vtr(
				MAC_Module* vtk,
                         	bool parallel ) 
//----------------------------------------------------------------------
{
   MAC_LABEL("FV_ParaviewPostProcessingWriter:: build_vtr") ;

   std::string pre = ( parallel ? "P" : "" ) ;
   
   vtk->add_entry( "type", MAC_String::create( vtk,
                                               pre+"RectilinearGrid" ) ) ;

   if( BINARY )
   {
     vtk->add_entry( "byte_order",
		      MAC_String::create( vtk, "LittleEndian" ) ) ;
     vtk->add_entry( "compressor",
      		MAC_String::create( vtk, "vtkZLibDataCompressor" ) ) ;
   }
   
   MAC_Module* grid = MAC_Module::create( vtk, pre+"RectilinearGrid" ) ;

   if( parallel )
     grid->add_entry( "WholeExtent", MAC_IntVector::
   					create( grid, WHOLE_EXTENT ) ) ;
   else
     grid->add_entry( "WholeExtent", MAC_IntVector::
   					create( grid, EXTENT ) ) ;
   
   MAC_Module* piece = MAC_Module::create( grid, "Piece" ) ;
   piece->add_entry( "Extent", MAC_IntVector::
   					create( piece, EXTENT ) ) ;
   MAC_Module* base = ( parallel ? grid : piece ) ;
   
   PRIMARY_GRID->write_grid(base, parallel) ;
   
   // Create PointData and CellData modules
   MAC_Module* PointData = MAC_Module::create( base, pre+"PointData" ) ;
   MAC_Module* CellData = MAC_Module::create( base, pre+"CellData" ) ;
   
   // Write fields in vtk  files
   list<FV_DiscreteField const*>::const_iterator il;
   for (il=FVFIELDS.begin(); il!=FVFIELDS.end() ; il++)
     (*il)->write_field( PointData, CellData );

   base->add_module( PointData ) ;
   base->add_module( CellData ) ;

   if( !parallel ) grid->add_module(piece) ;
   vtk->add_module( grid ) ;   
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: write_vtk( MAC_ModuleExplorer* vtk,
                           std::ofstream& file,
                           size_t level,
                           bool parallel )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: write_vtk" ) ;

   std::string data_array = ( parallel ? "PDataArray" : "DataArray" ) ;
   
   if( level == 0 )
      file << "<?xml version=\"1.0\"?>" << endl ;
   std::string bl( 2*level, ' ' ) ;
   bl = "\n" + bl ;
   std::string name =  vtk->name() ;
   if( name.find( "#" )<name.length() )
      name = name.substr( 0, name.find( "#" ) ) ;
   
   file << bl << "<" << name ;
   for( vtk->start_entry_iterator() ;
        vtk->is_valid_entry() ;
        vtk->go_next_entry() ) 
   {
      MAC_Data* data = vtk->data( 0 ) ;
      MAC_Data::Type dt = data->data_type() ;
      
      if( dt==MAC_Data::String  ) 
      {
         file << " " << vtk->keyword() << "=" ;
         data->print(file,0) ;
         file<<" " ;
      }
      else if( dt==MAC_Data::Int ) 
      {
         file << " " << vtk->keyword() << "=\"" ;
         data->print(file,0) ;
         file<<"\" " ;
      }
      else if( dt==MAC_Data::IntVector ) 
      {
         file << " " << vtk->keyword() << "=\"" ;
         intVector const& tab = data->to_int_vector() ;
	 for( size_t i=0 ; i<tab.size() ; i++ )
	   file << " " << tab(i) ;
        file<<"\" " ;
      }
      data->destroy() ;
   }
   file << ">" << endl;
   
   for( vtk->start_entry_iterator() ;
        vtk->is_valid_entry() ;
        vtk->go_next_entry() ) 
   {
      MAC_Data* data = vtk->data( 0 ) ;
      MAC_Data::Type dt = data->data_type() ;
      if( dt==MAC_Data::DoubleArray2D || dt==MAC_Data::DoubleVector ) 
      {
         file << "<" << data_array << " " ;        
         file << "Name=\"" << vtk->keyword() << "\" " ;

         std::string format = "ascii" ;
         
         if( BINARY && !parallel )
         {
            format = "appended" ;
            file << "offset=\"" << OFFSET << "\" " ;
         }
         
         if( dt == MAC_Data::DoubleArray2D )
         {
            int nbc = data->to_double_array2D().index_bound(0) ;
            if( nbc > 1 )
               file << "NumberOfComponents=\"3\" " ;
         }
         
         file << "type=\"Float32\" " ;
         file << "format=\"" << format << "\" " ;
         file << ">" << endl ;

         if( !parallel )
         {           
            if( dt==MAC_Data::DoubleArray2D ) 
            {
               doubleArray2D const& tab = data->to_double_array2D() ;
               size_t nbc = data->to_double_array2D().index_bound(0) ;
               size_t ncmax = ( nbc==1 ? 1 : 3 ) ;
            
	       start_output( sizeof_Float32, ncmax*tab.index_bound(1) ) ;
	       for( size_t i=0 ; i<tab.index_bound(1) ; i++ )
               {
                  for( size_t j=0 ; j<ncmax ; j++ )
		    write_double( file, ( j<nbc ? tab(j,i) : 0.0 ) ) ;
               }
	       flush(file) ;
               
            }
            else if( dt==MAC_Data::DoubleVector  )
            {
               doubleVector const& tab = data->to_double_vector() ;
	       if (tab.size() > 0)
	       {
	         start_output( sizeof_Float32, tab.size() ) ;
	         for( size_t i=0 ; i<tab.size() ; i++ )
	           write_double( file, tab(i) ) ;
		 flush(file) ;
	       }
            }
         }
         
         file << "</" << data_array << ">" << endl;
      }
      else if( ! ( dt==MAC_Data::String || dt==MAC_Data::Int
      				|| dt==MAC_Data::IntVector) ) 
      {
         MAC_Error::object()->raise_internal(
            " Bad type "+data->data_type() ) ;
      }
      data->destroy() ;
   }

   for( vtk->start_module_iterator() ;
        vtk->is_valid_module() ;
        vtk->go_next_module() ) 
   {
      MAC_ModuleExplorer* sexp = vtk->create_subexplorer( 0 ) ;
      write_vtk(sexp,file,level+1,parallel) ;
      sexp->destroy() ;
   }
   
   if( level == 0 && BINARY && !parallel ) 
   {
      file << bl << "<AppendedData encoding=\"raw\" > " << endl << "    _" ;
      file.write( BUFFER, OFFSET ) ;
      file << bl << "</AppendedData>" << endl ;
      delete [] BUFFER ; BUFFER = 0 ;
      ALLOCATED = 0 ;
      OFFSET = 0 ;
   }
   
   file << bl << "</" << name << ">" << endl;   
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: start_output( size_t size,
							size_t number )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: start_output" ) ;
   if( BINARY )
   {
      int current_output_size = size*number ;
      unsigned long ncomp = current_output_size 
      			+ (current_output_size+999)/1000 + 12 + sizeof_Int32 ;
      check_allocated( ncomp ) ;
      CURRENT_LENGTH = store_int(current_output_size) ;
   }
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: write_double( std::ofstream& file,
							double val )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: write_double" ) ;
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
void
FV_ParaviewPostProcessingWriter:: write_int( std::ofstream& file, int val )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: write_int" ) ;
   MAC_CHECK( sizeof_Int32==sizeof(int) ) ;
   
   if( BINARY )
      store_int(val) ;
   else
      file << " " << val ;
}




//-----------------------------------------------------------------------
size_t
FV_ParaviewPostProcessingWriter:: store_int( int val )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: store_int" ) ;
   MAC_CHECK( OFFSET + sizeof_Int32 <= ALLOCATED ) ;
   
   size_t result = OFFSET ;
   *((int*)&(BUFFER[OFFSET])) = val ;
   OFFSET += sizeof_Int32  ;
   return result ;
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: check_allocated( size_t size )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: check_allocated" ) ;
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
FV_ParaviewPostProcessingWriter:: flush( std::ofstream& file )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: flush" ) ;
  
   if( BINARY )
     compress_segment(CURRENT_LENGTH) ;
   file << endl ;
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: compress_segment( size_t seg )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: compress_segment" ) ;
   static size_t BlockSize = 32768 ;
   size_t size = (size_t)(*((int*)&BUFFER[seg])) ;
   
   size_t numFullBlocks = size / BlockSize;
   size_t lastBlockSize = size % BlockSize;
   size_t numBlocks = numFullBlocks + (lastBlockSize?1:0);

   size_t headerLength = numBlocks+3;

   int * CompressionHeader = new int[headerLength];
   CompressionHeader[0] = numBlocks;
   CompressionHeader[1] = BlockSize;
   CompressionHeader[2] = lastBlockSize;

   unsigned long encoded_buff_size = MAC::max(BlockSize,size)  ;
   unsigned char* encoded_buff = new unsigned char [ encoded_buff_size ] ;
   size_t encoded_offset = 0 ;
   for( size_t block=0 ; block<numBlocks ; block++ )
   {
      size_t buffer_start = seg + sizeof_Int32 + block*BlockSize ;
      size_t length = ( block+1<numBlocks || 
      	!lastBlockSize ? BlockSize : lastBlockSize ) ;
      unsigned char* to_encode = (unsigned char *)(&BUFFER[buffer_start]) ;
      unsigned char* encoded = &encoded_buff[encoded_offset] ;
      unsigned long ncomp = encoded_buff_size - encoded_offset ;
      
      if(compress2((Bytef*)encoded,
                   &ncomp,
                   (const Bytef*)to_encode,
                   length,
                   Z_DEFAULT_COMPRESSION) != Z_OK)
      {
         MAC_Error::object()->raise_plain(
            "Zlib error while compressing data.");
      }
      CompressionHeader[3+block] = ncomp ;
      encoded_offset += ncomp ;
   }
   
   OFFSET = seg ;
   check_allocated( headerLength*sizeof_Int32 + encoded_offset ) ;
   
   for(size_t i=0 ; i<headerLength ; i++ )
      store_int(CompressionHeader[i]) ;     

   for(size_t i=0 ; i<encoded_offset ; i++ )
      BUFFER[OFFSET++] = encoded_buff[i] ;

   if( OFFSET%4 != 0 )
      OFFSET = 4*( OFFSET/4 +1 ) ; // Re-allignement
   
   delete [] CompressionHeader ;
   delete [] encoded_buff ;
}




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: clearResultFiles( void )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: clearResultFiles" ) ;

   if ( MAC_Exec::communicator()->rank() == 0 )
   {
     FV::out() << "   CLEAR PARAVIEW RESULTS FILES IN DIRECTORY " << 
     	RES_DIRECTORY << endl;

     char* path = getenv("MAC_HOME");
     string exe(path);
     exe = exe+"/tools/ShellScripts/FVParaview_clear_exec "+RES_DIRECTORY;
     system(exe.c_str()); 
   }
}




//-----------------------------------------------------------------------
size_t
FV_ParaviewPostProcessingWriter:: getPreviousCycleNumber( void )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: getPreviousCycleNumber" ) ; 
   
   string tline,previous_tline,buffer,part;

   if ( MAC_Exec::communicator()->rank() == 0 )  
      FV::out() << "   GET PARAVIEW PREVIOUS CYCLE NUMBER ";   

   // Read pvd file
   ifstream fileIN( PVD_FILENAME.c_str(), ios::in );
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

   if ( MAC_Exec::communicator()->rank() == 0 )  
     FV::out() << CYCLE_NUMBER << endl;  
  
   return CYCLE_NUMBER ;     
}  




//-----------------------------------------------------------------------
void
FV_ParaviewPostProcessingWriter:: readTimeFile( FV_TimeIterator const* t_it,
      			size_t& cycle_number )  
//-----------------------------------------------------------------------
{
   MAC_LABEL( "FV_ParaviewPostProcessingWriter:: readTimeFile" ) ; 

   string tline, keyword ;
   double starting_time = t_it->time(), time_step = t_it->time_step(),
   	output_time ; 
   list<double> output_times ;
   list<string> output_time_lines ;	
   
   ifstream fileIN( PVD_FILENAME.c_str(), ios::in );
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
     PVD_STRINGS.append( *il );      
}   
