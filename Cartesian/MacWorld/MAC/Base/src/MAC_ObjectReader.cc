#include <MAC_ObjectReader.hh>

#include <MAC_Communicator.hh>
#include <MAC_Data.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_assertions.hh>
#include <MAC.hh>

#include <fstream>
#include <sstream>

using std::endl ;
using std::ifstream ;
using std::ostringstream ;
using std::istringstream ;
using std::string ;
using std::stack ;

struct MAC_ObjectReader_ERROR
{
   static void n0( string const& fname ) ;
   static void n1( string const& fname ) ;
   static void n2( size_t cycle_number, string const& name ) ;
   static void n3( string const& name ) ;
   static void n4( string const& class_name, string const& nn ) ;
   static void n5( void ) ;
   static void n6( string const& class_name ) ;
   static void n7( size_t stored_nb_rank, size_t nb_ranks ) ;
   static void n8( void ) ;   
} ;



//---------------------------------------------------------------------------
MAC_ObjectReader*
MAC_ObjectReader:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: create" ) ;

   MAC_ObjectReader* result = new MAC_ObjectReader( a_owner, exp ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( !result->positioned_in_a_valid_cycle() ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
MAC_ObjectReader:: MAC_ObjectReader( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , ROOT_MOD( 0 )
   , HEADER_MOD( 0 )
   , IFILE_NAME( exp->has_entry( "file_name" ) ?
	exp->string_data( "file_name" ) : MAC::undefined_string )
   , NB_CYCLES( 0 )
   , LAST_CYCLE( 0 )
   , iOBJECT( 0 )
   , initial_time( 0. )
{
   // Get restart file either from file name or from restart table file
   if ( exp->has_entry( "file_name" ) )
     IFILE_NAME = exp->string_data( "file_name" ) ;
   else if ( exp->has_entry( "restart_table_file_directory" ) )
   {
     string dir = exp->string_data( "restart_table_file_directory" );
     string rftable = dir + "/" + "RFTable.txt";

     ifstream fileIN( rftable.c_str(), std::ios::in  ) ;
     if( !fileIN ) MAC_ObjectReader_ERROR::n0( rftable ) ;               
     
     string tline;
     getline( fileIN, tline );
     while ( !fileIN.eof() ) 
     { 
       istringstream iss( tline );
       iss >> initial_time >> IFILE_NAME;
       getline( fileIN, tline );    
     }
     fileIN.close() ;
     IFILE_NAME = dir + "/" + IFILE_NAME ;      
   }
   else MAC_ObjectReader_ERROR::n8() ;

   // Get header module file name
   string dir = MAC::extract_root_directory( IFILE_NAME );
   string rootfilename = MAC::extract_root_file_name( IFILE_NAME );
   string twolast = rootfilename.substr( rootfilename.size() - 2 );
   if ( twolast == "#A" || twolast == "#B" )
     rootfilename = rootfilename.substr( 0, rootfilename.size() - 2 );   
   HEADERFILE = dir + "/" + rootfilename + "_header.mac";
   
   // Add rank in case of parallel
   string fname = IFILE_NAME ;
   MAC_Communicator const* com = MAC_Exec::communicator() ;
   if( com->nb_ranks()>1 )
   {
      ostringstream rank ;
      rank << "." << com->rank() ;
      fname  += rank.str() ;
   }   

   ifstream file( fname.c_str(), std::ios::in  ) ;
   if( !file ) MAC_ObjectReader_ERROR::n0( fname ) ;
   file.close() ;
   
   ifstream file2( HEADERFILE.c_str(), std::ios::in  ) ;
   if( !file2 ) MAC_ObjectReader_ERROR::n0( HEADERFILE ) ;
   file2.close() ;   

   ROOT_MOD = MAC_Module::create( this, "ROOT", fname ) ;
   HEADER_MOD = MAC_Module::create( this, "HEADER", HEADERFILE ) ;

   // Check fname structure :
   MAC_ModuleIterator* it = ROOT_MOD->create_module_iterator( 0 ) ;
   if( !it->is_valid() ) 
      MAC_ObjectReader_ERROR::n1( fname ) ;
   if( !( it->item()->name()=="communicator") ) 
      MAC_ObjectReader_ERROR::n1( fname ) ;
   for( it->go_next() ; it->is_valid() ; it->go_next() )
   {
      MAC_Module const* mm = it->item() ;
      if( !(mm->name().substr(0,6)=="cycle#") )
         MAC_ObjectReader_ERROR::n1( fname ) ;
      if( !mm->has_entry( "cycle_number" ) )
         MAC_ObjectReader_ERROR::n3( mm->name() ) ;
      int i = mm->data_of_entry( "cycle_number" )->to_int() ;
      if( i<0 ) MAC_ObjectReader_ERROR::n3( mm->name() ) ;
      LAST_CYCLE = (size_t) i ;
      ++NB_CYCLES ;
   }
   it->destroy() ;
   
   // Check HEADERFILE structure :
   it = HEADER_MOD->create_module_iterator( 0 ) ;
   if( !it->is_valid() ) 
      MAC_ObjectReader_ERROR::n1( HEADERFILE ) ;
   if( !( it->item()->name()=="header") ) 
      MAC_ObjectReader_ERROR::n1( HEADERFILE ) ;
   it->destroy() ;   

   // Check communicator :
   MAC_ModuleExplorer const* com_exp =
      MAC_ModuleExplorer::create( 0, ROOT_MOD->module( "communicator" ) ) ;
   size_t const nb_ranks = com_exp->int_data( "nb_ranks" ) ;
   size_t const rank = com_exp->int_data( "rank" ) ;
   if( nb_ranks!=com->nb_ranks() )
      MAC_ObjectReader_ERROR::n7( nb_ranks, com->nb_ranks() ) ;
   if( rank!=com->rank() )
      MAC_ObjectReader_ERROR::n1( fname ) ;
   com_exp->destroy() ; 
}




//---------------------------------------------------------------------------
MAC_ObjectReader:: ~MAC_ObjectReader( void )
//---------------------------------------------------------------------------
{
}




//---------------------------------------------------------------------------
MAC_Module*
MAC_ObjectReader:: header_module( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: header_module" ) ;

   // the existence of such a module has been tested in the constructor 
   // Not entirely clear why we need to create an iterator and take the 
   // first item, and not simply write MAC_Module* result = HEADER_MOD
   // as HEADER_MOD as a single main module called MODULE header
   MAC_ModuleIterator* it = HEADER_MOD->create_module_iterator( 0 ) ;
   MAC_Module* result = it->item() ;      	 
   it->destroy() ;      

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
size_t
MAC_ObjectReader:: nb_cycles( void ) const
//---------------------------------------------------------------------------
{
   size_t result = NB_CYCLES ;
   return( result ) ;
}




//---------------------------------------------------------------------------
std::string
MAC_ObjectReader:: input_file_name( void ) const
//---------------------------------------------------------------------------
{
   return( IFILE_NAME ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectReader:: seek_cycle( size_t cycle_number )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: seek_cycle" ) ;

   if( cycle_number==0 ) cycle_number=LAST_CYCLE ;

   ostringstream name ;
   name << "cycle#" << cycle_number ;
   if( !ROOT_MOD->has_module( name.str() ) ) 
      MAC_ObjectReader_ERROR::n2( cycle_number, name.str() ) ;

   MAC_Module* mm = ROOT_MOD->module( name.str() ) ;
   MODS.push( mm ) ;
   MOD_ITS.push( mm->create_module_iterator( this ) ) ; 

   if( mm->data_of_entry( "cycle_number" )->to_int()!=(int)cycle_number )
      MAC_ObjectReader_ERROR::n3( name.str() ) ;

   iOBJECT = 0 ;

   MAC_CHECK_POST( current_object_number() == 0 ) ;
   MAC_CHECK_POST( positioned_in_a_valid_cycle() ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ObjectReader:: positioned_in_a_valid_cycle( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: positioned_in_a_valid_cycle" ) ;

   bool result = !MOD_ITS.empty() ;
   return( result ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectReader:: close_cycle( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: close_cycle" ) ;
   MAC_CHECK_PRE( positioned_in_a_valid_cycle() ) ;

   iOBJECT = 0 ;

   MAC_ModuleIterator* it = MOD_ITS.top() ;
   MOD_ITS.pop() ;
   destroy_possession( it ) ;

   MODS.pop() ;

   if( !MODS.empty() ) MAC_ObjectReader_ERROR::n5() ;

   MAC_CHECK_POST( !positioned_in_a_valid_cycle() ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectReader:: start_object_retrieval( std::string const& class_name )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: start_object_retrieval" ) ;
   MAC_CHECK_PRE( positioned_in_a_valid_cycle() ) ;

   if( !MOD_ITS.top()->is_valid() ) MAC_ObjectReader_ERROR::n6( class_name ) ;

   MAC_Module* mm = MOD_ITS.top()->item() ;
   MODS.push( mm ) ;
   MOD_ITS.push( mm->create_module_iterator( this ) ) ; 

   MAC_ModuleExplorer const* exp = MAC_ModuleExplorer::create( 0, mm ) ;
   string const& nn = exp->string_data( "class" ) ;
   if( nn != class_name ) MAC_ObjectReader_ERROR::n4( class_name, nn ) ;

   iOBJECT = exp->int_data( "object_number" ) ;
   if( ! ( iOBJECT > 0 ) ) 
      MAC_Error::object()->raise_bad_data_value( exp, 
                                                 "object_number", 
                                                 "greater or equal to 1" ) ;
   exp->destroy() ;

   MAC_CHECK_POST( current_object_number() != 0 ) ;
}




//---------------------------------------------------------------------------
size_t
MAC_ObjectReader:: current_object_number( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: current_object_number" ) ;

   return( (size_t)iOBJECT ) ;
}




//---------------------------------------------------------------------------
double
MAC_ObjectReader:: get_initial_time( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: get_initial_time" ) ;

   return( initial_time ) ;
}




//---------------------------------------------------------------------------
bool
MAC_ObjectReader:: has_entry( std::string const& keyword ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: has_entry" ) ;
   MAC_CHECK_PRE( positioned_in_a_valid_cycle() ) ;
   MAC_CHECK_PRE( current_object_number() != 0 ) ;
   
   return( MODS.top()->has_entry( keyword ) ) ;
}




//---------------------------------------------------------------------------
MAC_Data const*
MAC_ObjectReader:: data_of_entry( std::string const& keyword ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: data_of_entry" ) ;
   MAC_CHECK_PRE( positioned_in_a_valid_cycle() ) ;
   MAC_CHECK_PRE( current_object_number() != 0 ) ;

   MAC_Data const* result = MODS.top()->data_of_entry( keyword ) ;

   MAC_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectReader:: end_object_retrieval( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: end_object_retrieval" ) ;

   MAC_ModuleIterator* it = MOD_ITS.top() ;
   MOD_ITS.pop() ;
   destroy_possession( it ) ;

   MODS.pop() ;

   MOD_ITS.top()->go_next() ;
}




//---------------------------------------------------------------------------
std::string
MAC_ObjectReader:: next_object_name_in_current_module( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: next_object_name" ) ;
   MAC_CHECK_PRE( positioned_in_a_valid_cycle() ) ;
   
   std::string result = MAC::undefined_string ;

   if ( MOD_ITS.top()->is_valid() )
   {
     MAC_Module* mm = MOD_ITS.top()->item() ;
     MAC_ModuleExplorer const* exp = MAC_ModuleExplorer::create( 0, mm ) ;
   
     if ( exp->has_entry( "name" ) )
       result = exp->string_data( "name" );
   
     exp->destroy() ;
   }   
   
   return( result ) ;
}




//---------------------------------------------------------------------------
std::string
MAC_ObjectReader:: next_object_class_in_current_module( void ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: next_object_name" ) ;
   MAC_CHECK_PRE( positioned_in_a_valid_cycle() ) ;
   
   std::string result = MAC::undefined_string ;

   if ( MOD_ITS.top()->is_valid() )
   {
     MAC_Module* mm = MOD_ITS.top()->item() ;
     MAC_ModuleExplorer const* exp = MAC_ModuleExplorer::create( 0, mm ) ;
   
     if ( exp->has_entry( "class" ) )
       result = exp->string_data( "class" );
   
     exp->destroy() ;
   }   
   
   return( result ) ;
}




//---------------------------------------------------------------------------
void
MAC_ObjectReader:: print( std::ostream& os, size_t indent_width ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ObjectReader:: print" ) ;

   std::string s( indent_width, ' ' ) ;
   os << s << "*** Initializing objects from reload files" << endl << endl;
   os << s << "      Header file name = " << HEADERFILE << endl;
   os << s << "      Object file name = " << IFILE_NAME << endl << endl;   
}




//---------------------------------------------------------------------------
bool
MAC_ObjectReader:: invariant( void ) const
//---------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( MODS.size() == MOD_ITS.size() ) ;
   return( true ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n0( string const& fname )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   unable to open file \"" << fname << "\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n1( string const& fname )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   file \"" << fname << "\" has an invalid structure" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n2( size_t cycle_number, string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   object retrieval from cycle " << cycle_number 
        << " is impossible " << endl 
        << "   since the underlying file has no module called \"" 
        << name << "\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n3( string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   module \"" << name << "\" has a missing or invalid" << endl
        << "   entry of keyword \"cycle_number\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n4( string const& class_name, string const& nn )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   attempt to retrieve an object of" << endl 
        << "   class \"" << class_name << "\" from a module whose " << endl
        << "   entry of keyword \"class_name\" is \"" << nn << "\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n5( void )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   impossible to close a cycle when an " << endl 
        << "   object retrieval is in progress" << endl
        << "   (\"start_object_retrieval\" and \"end_object_retrieval\"" 
        << endl
        << "    should be called the same number of times between two calls to"
	<< endl
	<< "    \"seek_cycle\" and \"close_cycle\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n6( string const& class_name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   attempt to retrieve an object of" << endl 
        << "   class \"" << class_name << "\" from a module who " << endl
        << "   does not contain any more object storage"  ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n7( size_t stored_nb_rank, size_t nb_ranks )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   object retrieval is impossible on "
        << nb_ranks << " processes" << endl 
        << "   because the storage has been done on "
        << stored_nb_rank << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//internal--------------------------------------------------------------
void 
MAC_ObjectReader_ERROR:: n8( void )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "MAC_ObjectReader :" << endl
        << "   missing one of the 2 keywords \"file_name\" or " 
        << "\"restart_table_file_directory\"" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
