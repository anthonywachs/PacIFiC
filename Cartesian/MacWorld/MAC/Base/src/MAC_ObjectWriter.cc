#include <MAC_ObjectWriter.hh>

#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_KeywordDataIterator.hh>
#include <MAC_Int.hh>
#include <MAC_Error.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_ListIterator.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_String.hh>
#include <MAC.hh>
#include <MAC_Application.hh>

#include <fstream>
#include <sstream>

struct MAC_ObjectWriter_ERROR
{
   static void n0( void ) ;
   static void n1( MAC_ModuleExplorer const* exp, std::string const& n ) ;
   static void n2( std::string const& n ) ;
   static void n3( std::string const& root0, std::string const& root1 ) ;   
} ;

//---------------------------------------------------------------------------
MAC_ObjectWriter*
MAC_ObjectWriter:: create( MAC_Object* a_owner,
                           MAC_ModuleExplorer const* exp,
                           MAC_ModuleExplorer const* header_exp )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   std::string const& w_type = exp->string_data( "type" ) ;
   MAC_ObjectWriter::MAC_ObjectWriterType writer_type =
                                        MAC_ObjectWriter::last_two_cycles ;
   
   if( w_type == "all_cycles_in_one_file" )
   {
      writer_type = MAC_ObjectWriter::all_cycles ;
   }
   else if( w_type == "cycles_in_separate_files" )
   {
      writer_type = MAC_ObjectWriter::per_one_cycle ;
   }
   else if( w_type == "last_two_cycles" )
   {
      writer_type = MAC_ObjectWriter::last_two_cycles ;
   }
   else
   {
      MAC_Error::object()->raise_bad_data_value(
         exp,
         "type",
         "   - \"all_cycles_in_one_file\"\n"
         "   - \"cycles_in_separate_files\"\n"
         "   - \"last_two_cycles\"" ) ;
   }
   
   MAC_ObjectWriter* result = new MAC_ObjectWriter( a_owner, writer_type,
                                                    exp, header_exp ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( !result->has_an_opened_cycle() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_ObjectWriter:: MAC_ObjectWriter( MAC_Object* a_owner,
                                     MAC_ObjectWriterType const writer_type,
                                     MAC_ModuleExplorer const* exp,
                                     MAC_ModuleExplorer const* header_exp )
//---------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , TYPE( writer_type )
   , OFILE_FORMAT( exp->has_entry( "output_format" ) ?
	exp->string_data( "output_format" ) : "hybrid" )
   , HEADER_EXP( header_exp != 0 ? header_exp->create_clone( this ) : 0 )
   , OFILE_NAME( )
   , OFILE_NAME0( )
   , OFILE_NAME1( )
   , RFTABLE( "RFTable.txt" )
   , HEADERFILE( )
   , B_WRITE_LAST_ITER( exp->has_entry( "force_write_at_last_time" ) ?
	exp->bool_data( "force_write_at_last_time" ) : false )
   , iCYCLE( 0 )
   , NB_OBJECTS( 0 )
{
   if( exp->has_entry( "output_format" ) )
   {
      if( OFILE_FORMAT!="text" && OFILE_FORMAT!="hybrid" )
      {
	 MAC_Error::object()->raise_bad_data_value(
	    exp,
	    "output_format",
	    "   - \"text\"\n   - \"hybrid\"" ) ;
      }
   }

   // The RFTable.txt keeps track of the restart files written versus time
   // RFTable.txt is written where restart files are written
   // Header file has _header.mac extension by construction, cannot be changed 
   // by users
   std::string RootSaving,filename,filename0,filename1,headerfilename;
   MAC_Communicator const* com = MAC_Exec::communicator() ;

   if ( TYPE == MAC_ObjectWriter::all_cycles )
   {
      OFILE_NAME = exp->string_data( "file_name" ) ;
      if( OFILE_NAME.empty() ) MAC_ObjectWriter_ERROR:: n1( exp, "file_name" ) ;
      
      size_t pos = OFILE_NAME.find_last_of( "/" );
      if ( pos == std::string::npos ) 
      {
        RootSaving = "."; 
	filename = OFILE_NAME;
      }
      else
      {
        RootSaving = OFILE_NAME.substr( 0, pos );
	filename = OFILE_NAME.substr( pos + 1 );   
      }
	   
      RFTABLE = RootSaving + "/" + RFTABLE ;
      
      pos = filename.find_first_of( "." ); 
      if ( pos == std::string::npos )      
	HEADERFILE = filename+"_header.mac";
      else
      {
	HEADERFILE = filename.substr( 0, pos );
	HEADERFILE.insert(pos,"_header.mac");	
      }

      headerfilename = HEADERFILE;
      HEADERFILE = RootSaving + "/" + HEADERFILE;       
      
      if( com->nb_ranks() > 1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILE_NAME += rank.str() ;
      }            
   }
   else if ( TYPE == MAC_ObjectWriter::per_one_cycle )
   {
      OFILE_NAME0 = exp->string_data( "files_basename" ) ;
      if( OFILE_NAME0.empty() ) MAC_ObjectWriter_ERROR:: n1( exp, 
      	"files_basename" ) ;
      
      size_t pos = OFILE_NAME0.find_last_of( "/" );
      if ( pos == std::string::npos ) 
      { 
        RootSaving = ".";
	filename0 = OFILE_NAME0;
      } 
      else 
      {
        RootSaving = OFILE_NAME0.substr( 0, pos );
	filename0 = OFILE_NAME0.substr( pos + 1 );   
      }
      
      RFTABLE = RootSaving + "/" + RFTABLE ;
      
      pos = filename0.find_first_of( "." ); 
      if ( pos == std::string::npos )      
	HEADERFILE = filename0+"_header.mac";
      else
      {
	HEADERFILE = filename0.substr( 0, pos );
	HEADERFILE.insert(pos,"_header.mac");	
      }

      headerfilename = HEADERFILE;
      HEADERFILE = RootSaving + "/" + HEADERFILE;         	
   }
   else if ( TYPE == MAC_ObjectWriter::last_two_cycles )
   {
      // Last_two_cycles forces the following:
      // 1) both cycles are in the same directory
      // 2) both cycles have the same root name and are distinguished by 
      // an additional "A" and "B"     
      OFILE_NAME0 = exp->string_data( "files_basename" ) ; 
      if( OFILE_NAME0.empty() ) MAC_ObjectWriter_ERROR:: n1( exp, 
      	"files_basename" ) ;
	
      std::string rootfilename ;
      size_t pos = OFILE_NAME0.find_last_of( "/" );
      if ( pos == std::string::npos ) 
      { 
        RootSaving = "."; 
	rootfilename = OFILE_NAME0;
      }	
      else 
      { 
        RootSaving = OFILE_NAME0.substr( 0, pos );
	rootfilename = OFILE_NAME0.substr( pos + 1 ); 
      }
      
      RFTABLE = RootSaving + "/" + RFTABLE ;      
      
      pos = rootfilename.find_first_of( "." ); 
      if ( pos == std::string::npos ) 
      {
        OFILE_NAME0 = rootfilename+"#A";
        OFILE_NAME1 = rootfilename+"#B";
	HEADERFILE = rootfilename+"_header.mac";
      }
      else
      {
        OFILE_NAME0 = rootfilename;
	OFILE_NAME0.insert(pos,"#A");
        OFILE_NAME1 = rootfilename;
	OFILE_NAME1.insert(pos,"#B");
	HEADERFILE = rootfilename.substr( 0, pos );
	HEADERFILE.insert(pos,"_header.mac");	
      }
      
      filename0 = OFILE_NAME0;
      filename1 = OFILE_NAME1;
      headerfilename = HEADERFILE;
      
      OFILE_NAME0 = RootSaving + "/" + OFILE_NAME0;
      OFILE_NAME1 = RootSaving + "/" + OFILE_NAME1; 
      HEADERFILE = RootSaving + "/" + HEADERFILE;     
            
      if( com->nb_ranks() > 1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILE_NAME0 += rank.str() ;
         OFILE_NAME1 += rank.str() ;
      }
      OFILE_NAME = OFILE_NAME1 ;
      // Comment: based on the above, the 1st file name, if not swapped,
      // will be OFILE_NAME0, see set_file_name()             	           
   }


   if ( com->rank() == 0 )   
   {
      MAC::out() << std::endl << "*** Saving files features for reload" 
      		<< std::endl << std::endl;
      MAC::out() << "   Output format = " << OFILE_FORMAT << std::endl;      
      MAC::out() << "   Directory = " << RootSaving << std::endl;
      MAC::out() << "   Type = ";
      if ( TYPE == MAC_ObjectWriter::all_cycles )
      {
        MAC::out() << "all_cycles_in_one_file" << std::endl;
        MAC::out() << "   file_name = " << filename << std::endl;      
      }	
      else if ( TYPE == MAC_ObjectWriter::per_one_cycle )
      {
        MAC::out() << "cycles_in_separate_files" << std::endl;
	MAC::out() << "   file_name = " << filename0 << std::endl;
      }
      else if ( TYPE == MAC_ObjectWriter::last_two_cycles )
      {      
        MAC::out() << "last_two_cycles" << std::endl;
	MAC::out() << "   file_name_0 = " << filename0 << std::endl;
        MAC::out() << "   file_name_1 = " << filename1 << std::endl;        	
      }
      MAC::out() << "   header = " << headerfilename << std::endl;
      MAC::out() << "   Table full file name = " << RFTABLE << std::endl;
      MAC::out() << std::endl << std::endl;	
   }   
}




//---------------------------------------------------------------------------
MAC_ObjectWriter:: ~MAC_ObjectWriter( void )
//---------------------------------------------------------------------------
{
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: start_cycle( double const& time )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: start_cycle" ) ;
   MAC_CHECK_PRE( !has_an_opened_cycle() ) ;
   MAC_SAVEOLD( size_t, cycle_number, cycle_number() ) ;
   
   ++iCYCLE ;
   
   initialize_saving_file() ;
   
   write_rftable_file( time ) ;

   std::ostringstream mn ;
   mn << "cycle#" << iCYCLE ;
   MAC_Module* mod = MAC_Module::create( this, mn.str() ) ;

   mod->add_entry( "cycle_number", MAC_Int::create( mod, iCYCLE ) ) ;

   MODS.push( mod ) ;

   MAC_CHECK_POST( has_an_opened_cycle() ) ;
   MAC_CHECK_POST( cycle_number() == OLD( cycle_number ) + 1 ) ;
   MAC_CHECK_POST( current_object_number() == 0 ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: terminate_cycle( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: terminate_cycle" ) ;
   MAC_CHECK_PRE( has_an_opened_cycle() ) ;

   MAC_Module* mod = MODS.top() ;
   MODS.pop() ;

   if( !MODS.empty() ) MAC_ObjectWriter_ERROR::n0() ;

   mod->write( OFILE_NAME, OFILE_FORMAT ) ;

   destroy_possession( mod ) ;

   NB_OBJECTS=0 ;

   MAC_CHECK_POST( !has_an_opened_cycle() ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ObjectWriter:: has_an_opened_cycle( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: has_an_opened_cycle" ) ;
   return( !MODS.empty() ) ;
}




//---------------------------------------------------------------------------
size_t
MAC_ObjectWriter:: cycle_number( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: cycle_number" ) ;
   return( iCYCLE ) ;
}




//----------------------------------------------------------------------
std::string
MAC_ObjectWriter:: get_current_file_name( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: get_current_file_name" ) ;

   return( OFILE_NAME );
}




//----------------------------------------------------------------------
bool
MAC_ObjectWriter:: force_write_at_last_time( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Application:: force_write_at_last_time" ) ;

   return( B_WRITE_LAST_ITER );
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: start_new_object( std::string const& class_name )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: start_new_object" ) ;
   MAC_CHECK_PRE( has_an_opened_cycle() ) ;
   MAC_SAVEOLD( size_t, current_object_number, current_object_number() ) ;

   ++NB_OBJECTS ;

   std::ostringstream name ;
   name << "object#" << NB_OBJECTS ;

   MAC_Module* mod = MAC_Module::create( MODS.top(), name.str() ) ;
   mod->add_entry( "class", MAC_String::create( mod, class_name ) ) ;
   mod->add_entry( "object_number", MAC_Int::create( mod, (int)NB_OBJECTS ) ) ;
   MODS.top()->add_module( mod ) ;
   MODS.push( mod ) ;

   MAC_CHECK_POST( current_object_number()==OLD(current_object_number) + 1 ) ;
}




//---------------------------------------------------------------------------
size_t
MAC_ObjectWriter:: current_object_number( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: current_object_number" ) ;
   MAC_CHECK_PRE( has_an_opened_cycle() ) ;

   return( MODS.size()-1 ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: add_entry( std::string const& keyword, MAC_Data* data )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: add_entry" ) ;
   MAC_CHECK_PRE( has_an_opened_cycle() ) ;
   MAC_CHECK_PRE( current_object_number() != 0 ) ;
   MAC_CHECK_PRE( data->owner() == 0 ) ;

   data->set_owner( MODS.top() ) ;
   MODS.top()->add_entry( keyword, data ) ;

   MAC_CHECK_POST( data->is_under_ownership_of( this ) ) ;
}




//---------------------------------------------------------------------------
void MAC_ObjectWriter:: finalize_object( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "void MAC_ObjectWriter:: finalize_object" ) ;
   MAC_CHECK_PRE( has_an_opened_cycle() ) ;
   MAC_SAVEOLD( size_t, current_object_number, current_object_number() ) ;

   MODS.pop() ;

   MAC_CHECK_POST( current_object_number() == OLD(current_object_number)-1 ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: initialize_saving_file( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: initialize_saving_file" ) ;

   static bool first = true ;
   
   if( first || TYPE != all_cycles )
   {
      set_file_name() ;
      
      std::ofstream file( OFILE_NAME.c_str(), std::ios::out | std::ios::trunc );
      if( !file ) MAC_ObjectWriter_ERROR:: n2( OFILE_NAME ) ;
      file.close() ;
      if( OFILE_FORMAT=="hybrid" )
      {
         std::string const bin_file_name =  OFILE_NAME+ ".bin" ;
         std::ofstream file_bin( bin_file_name.c_str(),
                                 std::ios::out |
                                 std::ios::binary | 
                                 std::ios::trunc ) ;
         if( !file_bin ) MAC_ObjectWriter_ERROR:: n2( bin_file_name ) ;
         file_bin.close() ;
      }

      std::ofstream os( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
      os.close() ;
      write_communicator() ;
   }
   if ( first ) write_header() ;
   first = false ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: set_file_name( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: set_file_name" ) ;

   if( TYPE == all_cycles )
   {
      // Nothing to do
   }
   else if( TYPE == per_one_cycle )
   {
      MAC_ASSERT( iCYCLE<99999 ) ;
      std::ostringstream tmp ;
      tmp << iCYCLE ;
      std::string nb_string = tmp.str(), rootname, ext ;
      size_t pos = OFILE_NAME0.find_first_of( "." ); 
      if ( pos == std::string::npos ) 
      {     
	rootname = OFILE_NAME0;
	ext = "";
      }	
      else
      {
	rootname = OFILE_NAME0.substr( 0, pos );
	ext = OFILE_NAME0.substr( pos );
      }      
      OFILE_NAME = rootname+".00000";
      OFILE_NAME.replace( OFILE_NAME.length()-nb_string.length(),
                          nb_string.length(), nb_string ) ;
      OFILE_NAME += ext ;
      MAC_Communicator const* com = MAC_Exec::communicator() ;
      if( com->nb_ranks()>1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILE_NAME += rank.str() ;
      }
   }
   else if( TYPE == last_two_cycles )
   {
      if( OFILE_NAME == OFILE_NAME0 )
      {
         OFILE_NAME = OFILE_NAME1 ;
      }
      else if( OFILE_NAME == OFILE_NAME1 )
      {
         OFILE_NAME = OFILE_NAME0 ;
      }
      // Comment: based on the above, the 1st file name, if not swapped,
      // will be OFILE_NAME0, from the constructor
   }
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: write_header( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: write_header" ) ;

   MAC_Communicator const* com = MAC_Exec::communicator() ;

   if ( com->rank() == 0 )
   {
      std::ofstream os( HEADERFILE.c_str(), std::ios::out | std::ios::trunc ) ;
      os << "MODULE header" << std::endl ;
      os.close() ;

      // The header is always written in text format
      if( HEADER_EXP != 0 )  HEADER_EXP->write( HEADERFILE, "text" ) ;

      os.open( HEADERFILE.c_str(), std::ios::out | std::ios::app ) ;
      os << "END MODULE header" << std::endl ;
      os.close() ;
   }
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: write_communicator( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: write_communicator" ) ;

   MAC_Communicator const* com = MAC_Exec::communicator() ;
   std::ofstream os( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
   os << "MODULE communicator" << std::endl ;
   os << "   nb_ranks = " << com->nb_ranks() << std::endl ;
   os << "   rank = " << com->rank() << std::endl ;
   os << "END MODULE communicator" << std::endl ;
   os.close() ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: write_rftable_file( double const& time ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: write_rftable_file" ) ;

   static bool first = true ;
   
   MAC_Communicator const* com = MAC_Exec::communicator() ;
   if ( com->rank() == 0 )
   {     
     std::ofstream os;
     if ( first && !MAC_Application::is_follow() )
       os.open( RFTABLE.c_str(), std::ios::out ) ;
     else
       os.open( RFTABLE.c_str(), std::ios::out | std::ios::app ) ;  
   
     std::string fname;
     size_t pos = OFILE_NAME.find_last_of( "/" );
     if ( pos == std::string::npos ) fname = OFILE_NAME; 
     else fname = OFILE_NAME.substr( pos + 1 );
     
     // Subtract the rank number at the end before writing in RFTable file
     if ( com->nb_ranks() > 1 )
     {
       pos = fname.find_last_of( "." );
       if ( pos != std::string::npos ) fname = fname.substr( 0, pos );       
     } 
   
     os << time << " " << fname << std::endl;
     os.close() ;
   }
   
   first = false ;   
}




//---------------------------------------------------------------------------
void
MAC_ObjectWriter:: swap_restart_file_names( std::string const& inputfilename ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectWriter:: swap_restart_file_names" ) ;

   if ( TYPE == MAC_ObjectWriter::last_two_cycles )
   {       
      MAC_Communicator const* com = MAC_Exec::communicator() ; 
      std::string opfname = OFILE_NAME0;
      size_t pos = 0;       

      if( com->nb_ranks() > 1 )
      {
        pos = OFILE_NAME0.find_last_of( "." );
	if ( pos != std::string::npos ) opfname = OFILE_NAME0.substr( 0, 
		pos );
      }
      
      pos = opfname.find_last_of( "/" );  
      if ( pos != std::string::npos ) opfname = 
      	opfname.substr( pos + 1 );
	
      std::string ipfname = inputfilename;
      pos = inputfilename.find_last_of( "/" );  
      if ( pos != std::string::npos ) ipfname = 
      	inputfilename.substr( pos + 1 );
      
      if ( opfname == ipfname ) OFILE_NAME = OFILE_NAME0;
      // Comment: based on the above, the 1st file name
      // will now be OFILE_NAME1, which is different from inputfilename   	
   }     
}



//---------------------------------------------------------------------------
bool
MAC_ObjectWriter:: invariant( void ) const
//---------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( MODS.size() <= NB_OBJECTS ) ;
   return( true ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectWriter_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "MAC_ObjectWriter :" << std::endl
        << "   impossible to terminate a cycle when an " << std::endl 
        << "   object storage is in progress" << std::endl
        << "   (\"start_new_object\" and \"finalize_object\"" 
        << std::endl
        << "    should be called the same number of times between two calls to"
	<< std::endl
	<< "    \"start_cycle\" and \"terminate_cycle\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectWriter_ERROR:: n1( MAC_ModuleExplorer const* exp,
                             std::string const& n )
//internal--------------------------------------------------------------
{
   MAC_Error::object()->raise_bad_data_value(
      exp, n, "a no empty string is expected" ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectWriter_ERROR:: n2( std::string const& n )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "MAC_ObjectWriter :\n" ;
   mess += "   Saving failure : unable to open file \"" ;
   mess += n ;
   mess += "\" for writing" ;
   MAC_Error::object()->raise_plain( mess ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectWriter_ERROR:: n3( std::string const& root0, 
	std::string const& root1 )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "MAC_ObjectWriter :\n" ;
   mess += "   Last_two_cycles mode: saving directory \"" ;
   mess += root0 ;
   mess += "\" and saving directory \"" ;
   mess += root1 ;
   mess += "\" do not match.\n" ;
   mess += "Matching saving directories in path for file_name_0 and" ;
   mess += " file_name_1 is mandatory." ;  
   MAC_Error::object()->raise_plain( mess ) ;
}
